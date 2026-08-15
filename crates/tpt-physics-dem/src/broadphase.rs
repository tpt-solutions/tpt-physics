//! Uniform-grid spatial hashing for broad-phase neighbour queries.

use crate::particle::Particle;

/// A uniform cubic grid broad-phase accelerator.
///
/// Particles are binned into cells of side `cell_size` (typically twice the
/// largest particle radius); two particles can only collide if their centres lie
/// in the same or adjacent cells.
pub struct SpatialHash {
    #[allow(dead_code)]
    cell_size: f64,
    /// `cell_key -> list of particle indices`.
    cells: std::collections::HashMap<[i64; 3], Vec<usize>>,
}

impl SpatialHash {
    /// Build a hash from `particles` using the given `cell_size`.
    pub fn build(particles: &[Particle], cell_size: f64) -> Self {
        let mut cells: std::collections::HashMap<[i64; 3], Vec<usize>> = std::collections::HashMap::new();
        for (i, p) in particles.iter().enumerate() {
            let key = Self::key(p.position, cell_size);
            cells.entry(key).or_default().push(i);
        }
        SpatialHash { cell_size, cells }
    }

    fn key(p: [f64; 3], cell_size: f64) -> [i64; 3] {
        let f = 1.0 / cell_size;
        [
            (p[0] * f).floor() as i64,
            (p[1] * f).floor() as i64,
            (p[2] * f).floor() as i64,
        ]
    }

    /// Candidate collision pairs `(i, j)` with `i < j`, within the 27-cell
    /// neighbourhood of each particle.
    pub fn candidate_pairs(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for (key, idxs) in &self.cells {
            for dx in -1..=1i64 {
                for dy in -1..=1i64 {
                    for dz in -1..=1i64 {
                        let nk = [key[0] + dx, key[1] + dy, key[2] + dz];
                        if let Some(neigh) = self.cells.get(&nk) {
                            for &i in idxs {
                                for &j in neigh {
                                    if i < j {
                                        pairs.push((i, j));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_nearby_pairs_only() {
        // Three particles: 0 and 1 close, 2 far away.
        let ps = vec![
            Particle::new([0.0, 0.0, 0.0], [0.0; 3], 0.5, 1000.0),
            Particle::new([0.5, 0.0, 0.0], [0.0; 3], 0.5, 1000.0),
            Particle::new([100.0, 100.0, 100.0], [0.0; 3], 0.5, 1000.0),
        ];
        let h = SpatialHash::build(&ps, 1.0);
        let pairs = h.candidate_pairs();
        assert!(pairs.contains(&(0, 1)));
        assert!(!pairs.contains(&(0, 2)));
        assert!(!pairs.contains(&(1, 2)));
    }
}
