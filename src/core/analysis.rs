pub mod hclust {

    pub fn haversine(
        (lat1, lon1): (f64, f64),
        (lat2, lon2): (f64, f64),
    ) -> f64 {
        const EARTH_RADIUS: f64 = 3958.756; // miles

        let (lat1, lon1) = (lat1.to_radians(), lon1.to_radians());
        let (lat2, lon2) = (lat2.to_radians(), lon2.to_radians());

        let delta_lat = lat2 - lat1;
        let delta_lon = lon2 - lon1;
        let x = (delta_lat / 2.0).sin().powi(2)
            + lat1.cos() * lat2.cos() * (delta_lon / 2.0).sin().powi(2);
        2.0 * EARTH_RADIUS * x.sqrt().atan()
    }
}

#[cfg(test)]
mod tests {
    use bio::io::fasta::Record;
    use kodama::{linkage, Method};

    use rayon::iter::IntoParallelRefIterator;
    use crate::core::analysis::hclust::haversine;
    use crate::core::io::read_fasta;
    use crate::core::utils::{fasta_distance_jukes_cantor_number, remove_empty};
    use crate::core::utils::pairs::Pairwise;

    #[test]
    fn test_haversine() {
        // From our data set. Each coordinate pair corresponds to a single observation.
        let coordinates = [
            (42.5833333, -71.8027778),
            (42.2791667, -71.4166667),
            (42.3458333, -71.5527778),
            (42.1513889, -71.6500000),
            (42.3055556, -71.5250000),
            (42.2694444, -71.6166667),
        ];
        // Build our condensed matrix by computing the dissimilarity between all
        // possible coordinate pairs.
        let mut condensed = vec![];
        for row in 0..coordinates.len() - 1 {
            for col in row + 1..coordinates.len() {
                condensed.push(haversine(coordinates[row], coordinates[col]));
            }
        }
        dbg!(&condensed);

        // The length of a condensed dissimilarity matrix is always equal to
        // `N-choose-2`, where `N` is the number of observations.
        assert_eq!(
            condensed.len(),
            (coordinates.len() * (coordinates.len() - 1)) / 2
        );
        use kodama::{linkage, Method};

        let dend = linkage(&mut condensed, coordinates.len(), Method::Average);
        dbg!(&dend);
    }

    #[test]
    fn test_distance_clust() {
        use rayon::iter::*;
        let file = "resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa";
        match read_fasta(file) {
            Ok(x) => {
                let mut x: Vec<Record> =
                    remove_empty(x).iter().take(10).cloned().collect();
                let distances = x.pairwise_map_condensed_upper(
                    fasta_distance_jukes_cantor_number,
                );
                let mut values: Vec<f64> =
                    distances.par_iter().map(|x| x.x).collect();
                let dend = linkage(&mut values, 10, Method::Single);
                dbg!(&dend);
            }
            Err(e) => {
                panic!("Error: {}", e)
            }
        };
    }
}
