//! Real GPU compute path (feature `gpu`) using [`wgpu`].
//!
//! Implements the sparse matrix–vector product `y = A x` as a WGSL compute
//! shader dispatched on the first available GPU adapter. The operator is
//! converted to a dense `f32` matrix on the host (the GPU path targets large
//! problems; `f32` keeps the shader portable across Vulkan/Metal/DX12
//! adapters). When no adapter is found (headless CI, no Vulkan/Metal/DX12) the
//! routine returns [`SolverError::BackendUnavailable`] — the same graceful
//! fallback the CPU path historically returned — so callers degrade cleanly.
//!
//! This module is only compiled when the `gpu` feature is enabled:
//!
//! ```text
//! cargo build -p tpt-physics-solver --features gpu
//! ```

use crate::error::SolverError;
use pollster::FutureExt;
use tpt_fem_sparse::Csr;

/// WGSL compute shader: `y[i] = Σ_j A[i, j] · x[j]` over a dense `n×n` matrix
/// stored row-major in a storage buffer.
const SHADER: &str = r#"
@group(0) @binding(0) var<storage, read> a: array<f32>;
@group(0) @binding(1) var<storage, read> x: array<f32>;
@group(0) @binding(2) var<storage, read_write> y: array<f32>;
@group(0) @binding(3) var<uniform> n: u32;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= n) { return; }
    var s: f32 = 0.0;
    for (var j: u32 = 0u; j < n; j = j + 1u) {
        s = s + a[i * n + j] * x[j];
    }
    y[i] = s;
}
"#;

/// Convert a slice of `f32` into little-endian bytes without `unsafe`.
fn f32_to_le_bytes(values: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for &v in values {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

/// Perform `y = A x` on the GPU.
///
/// Returns [`SolverError::BackendUnavailable`] if no GPU adapter can be
/// acquired in this environment.
pub fn matvec_gpu(a: &Csr, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
    let n = a.nrows;
    if a.ncols != n {
        return Err(SolverError::NotSquare {
            nrows: a.nrows,
            ncols: a.ncols,
        });
    }
    if x.len() != n || y.len() != n {
        return Err(SolverError::NotSquare {
            nrows: x.len(),
            ncols: n,
        });
    }

    // Host-side CSR → dense f32 conversion.
    let mut dense = vec![0.0f32; n * n];
    for r in 0..n {
        for idx in a.row_ptrs[r]..a.row_ptrs[r + 1] {
            let c = a.col_ind[idx];
            dense[r * n + c] = a.values[idx] as f32;
        }
    }
    let xf: Vec<f32> = x.iter().map(|v| *v as f32).collect();

    let yf = compute(&dense, &xf, n).map_err(SolverError::BackendUnavailable)?;
    for i in 0..n {
        y[i] = yf[i] as f64;
    }
    Ok(())
}

/// Build a wgpu device, run the matvec kernel, and read back `y`.
fn compute(a: &[f32], x: &[f32], n: usize) -> Result<Vec<f32>, String> {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::all(),
        ..Default::default()
    });
    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
        })
        .block_on()
        .ok_or_else(|| {
            "no GPU adapter available (run on a machine with Vulkan/Metal/DX12, \
             or disable the `gpu` feature)"
                .to_string()
        })?;

    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor::default(), None)
        .block_on()
        .map_err(|e| format!("failed to acquire GPU device: {e:?}"))?;

    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("tpt-physics matvec"),
        source: wgpu::ShaderSource::Wgsl(SHADER.into()),
    });

    let bind_group_layout =
        device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("matvec bg"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

    let pipeline_layout =
        device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("matvec pl"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("matvec pipeline"),
        layout: Some(&pipeline_layout),
        module: &shader,
        entry_point: "main",
        compilation_options: Default::default(),
    });

    let a_bytes = f32_to_le_bytes(a);
    let x_bytes = f32_to_le_bytes(x);
    let n_bytes = (n as u32).to_le_bytes();

    let a_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("a"),
        size: a_bytes.len() as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let x_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("x"),
        size: x_bytes.len() as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let y_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("y"),
        size: (n * std::mem::size_of::<f32>()) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let n_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("n"),
        size: n_bytes.len() as u64,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    queue.write_buffer(&a_buf, 0, &a_bytes);
    queue.write_buffer(&x_buf, 0, &x_bytes);
    queue.write_buffer(&n_buf, 0, &n_bytes);

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("matvec"),
        layout: &bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: a_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: x_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: y_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: n_buf.as_entire_binding() },
        ],
    });

    let mut encoder =
        device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("matvec enc") });
    {
        let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("matvec pass"),
            timestamp_writes: None,
        });
        cpass.set_pipeline(&pipeline);
        cpass.set_bind_group(0, &bind_group, &[]);
        cpass.dispatch_workgroups((n as u32 + 63) / 64, 1, 1);
    }
    // Copy the result buffer to a mappable staging buffer.
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("y staging"),
        size: (n * std::mem::size_of::<f32>()) as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });
    encoder.copy_buffer_to_buffer(&y_buf, 0, &staging, 0, staging.size());
    queue.submit(Some(encoder.finish()));

    let slice = staging.slice(..);
    slice.map_async(wgpu::MapMode::Read, |_| {});
    // Wait for the map to complete (wgpu 0.20 `map_async` does not itself
    // return a future to drive).
    device.poll(wgpu::Maintain::Wait);
    // Reconstruct f32s from the little-endian bytes (safe: we wrote f32 LE).
    let view = slice.get_mapped_range();
    let mut result = vec![0.0f32; n];
    for (i, chunk) in view.chunks_exact(4).enumerate() {
        let mut b = [0u8; 4];
        b.copy_from_slice(chunk);
        result[i] = f32::from_le_bytes(b);
    }
    drop(view);
    staging.unmap();
    Ok(result)
}
