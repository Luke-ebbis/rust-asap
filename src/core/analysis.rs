/// Object for pairwise distance analysis
pub mod distance {
    use bio::io::fasta::Record;
    use derivative::Derivative;
    use num_traits::Signed;
    use crate::core::utils::{fasta_distance_jukes_cantor_number, remove_empty};
    use crate::core::utils::pairs::{Pair, Pairwise};
    use derive_builder::Builder;
    use kodama::Float;
    use crate::core::io::read_fasta;

    /// For now it only works with acutal (symetric) distances.
    #[derive(Builder)]
    struct DistanceAnalysis <T: num_traits::Signed + ?Sized + Clone + Send,
        F: Fn(Record, Record) -> Pair<Record, T>
        + Send
        + 'static
        + Sync > {

        // Mandatory Field:
        /// Data to calculate with.
        data: Vec<Record>,

        // Or you can set the default
        /// Distance function
        f: F,
    }


    impl<T: Signed + ?Sized + Clone + Send,
        F: Fn(Record, Record) -> Pair<Record, T>
        +  'static + Sync + std::marker::Send> DistanceAnalysis< T, F> where Vec<Pair<bio::io::fasta::Record, f64>>: From<Vec<Pair<bio::io::fasta::Record, T>>> {
        fn run(mut self) -> CondensedDistanceMatrix {
            CondensedDistanceMatrix {matrix: self.data.pairwise_map_condensed_upper(self.f).into()}
        }
    }

    // impl<T: Signed + ?Sized + Clone + Send,
    //      F: Fn(Record, Record) -> Pair<Record, T>
    //       +  'static + Sync + std::marker::Send> DistanceAnalysisBuilder <T, F> {
    //
    //     // Private helper method with access to the builder struct.
    //     fn default_f(&self) -> Result<F, String> {
    //         Ok(fasta_distance_jukes_cantor_number)
    //     }
    // }

    #[test]
    fn test_distance() {
        let recs = read_fasta("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa").unwrap();
        let recs = remove_empty(recs);
        let distanceAnalysis = DistanceAnalysisBuilder::create_empty().data(recs).f(fasta_distance_jukes_cantor_number).build().unwrap();
        let mat = distanceAnalysis.run();
        dbg!(&mat);
        let max = mat.max();
        dbg!(&max);
    }

    #[derive(Debug)]
    #[derive(PartialEq, PartialOrd, Clone)]
    struct CondensedDistanceMatrix {
        matrix: Vec<Pair<Record, f64>>,
    }

    impl CondensedDistanceMatrix
    {
        pub fn iter(&self) -> impl Iterator + use<'_> {
            // Logic to return an iterator over the internal collection
            self.matrix.iter()
        }
    }


    impl CondensedDistanceMatrix {

        fn remove_inf(self) -> Vec<Pair<Record, f64>> {
            self.matrix.into_iter().filter(|x: &Pair<Record, f64>| x.x != f64::INFINITY).collect()
        }

        /// Find the max value of the matrix that is not infinite
        pub fn max(self) -> Pair<Record, f64> {
            self.remove_inf().iter().max_by(|x,y| x.partial_cmp(y).unwrap()).cloned().unwrap()
        }

        /// Find the min value of the matrix that is not infinite
        pub fn min(self) -> Pair<Record, f64> {
            self.remove_inf().iter().min_by(|x,y| x.partial_cmp(y).unwrap()).cloned().unwrap()
        }
    }


    #[derive(Debug)]
    pub struct SequenceDistanceAnalysis<T: num_traits::Signed + ?Sized + Clone + Send,
        F: Fn(Record, Record) -> Pair<Record, T>
        + Send
        + 'static
        + Sync,
    > {
        distance: F
    }



    impl<T: Float + num_traits::Signed + std::marker::Send,
        F: Fn(Record, Record) -> Pair<Record, T>
        + Send
        + 'static
        + Sync,
    > SequenceDistanceAnalysis<T, F> where Vec<Pair<bio::io::fasta::Record, f64>>: From<Vec<Pair<bio::io::fasta::Record, T>>> {
        /// Construct a sequence distance analysis.
        ///
        /// ```
        /// use asap::core::analysis::distance::SequenceDistanceAnalysis;
        /// use asap::core::utils::fasta_distance_jukes_cantor_number;
        /// let sequence_analysis = SequenceDistanceAnalysis::new(fasta_distance_jukes_cantor_number);
        /// ```
        pub fn new(distance: F)  -> Self {
            SequenceDistanceAnalysis {distance }
        }


        /// Run using the specified distance function.
        pub fn run(self, records: Vec<Record>) -> CondensedDistanceMatrix {
            let mut fasta_distances = records.clone()
                .pairwise_map_condensed_upper(self.distance);
           CondensedDistanceMatrix { matrix: fasta_distances.into()}
        }
    }


}

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


    pub fn hclust() -> () {

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
