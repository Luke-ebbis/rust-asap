use bio::alignment::distance::levenshtein;
use bio::io::fasta::Record;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefIterator;
use crate::core::utils::distances::jukes_cantor;
use crate::core::utils::pairs::Pair;

/// A struct to deal with paired operations
pub mod pairs {
    use bio::io::fasta::Record;

    use num_traits::Signed;
    use rayon::iter::IntoParallelRefIterator;
    use rayon::iter::ParallelIterator;
    use std::cmp::Ordering;
    use std::collections::HashSet;

    // Condensed Dissimilarity matrix.
    // This will be needed for the Dendrogram.
    // pub struct PairwiseMatrix<T> {
    //     x: String,
    //     x: T,
    // }

    /// # A Pair struct
    ///
    /// An object to store two objects with a value.
    #[derive(Eq, Ord, Debug, Copy, Clone)]
    pub struct Pair<A, T: Signed + Clone>
    where
        A: Clone,
        T: Clone,
    {
        pub(crate) a: A,
        pub(crate) b: A,
        pub(crate) x: T,
    }

    impl<A: Clone, T: Signed + Clone> Pair<A, T> {
        pub fn new(
            a: A,
            b: A,
            x: T,
        ) -> Self {
            Pair { a, b, x }
        }
    }

    /// For a pair, the sorting happens on the `x` field.
    impl<A, T> PartialEq<Self> for Pair<A, T>
    where
        T: Signed + PartialEq + Clone,
        A: Clone,
    {
        fn eq(
            &self,
            other: &Self,
        ) -> bool {
            self.x == other.x
        }
    }
    impl<A, T> PartialOrd for Pair<A, T>
    where
        T: Signed + PartialEq + PartialOrd + Clone,
        A: Clone,
    {
        fn partial_cmp(
            &self,
            other: &Self,
        ) -> Option<Ordering> {
            self.x.partial_cmp(&other.x)
        }

        fn lt(
            &self,
            other: &Self,
        ) -> bool {
            self.x.lt(&other.x)
        }

        fn le(
            &self,
            other: &Self,
        ) -> bool {
            self.x.le(&other.x)
        }

        fn gt(
            &self,
            other: &Self,
        ) -> bool {
            self.x.gt(&other.x)
        }

        fn ge(
            &self,
            other: &Self,
        ) -> bool {
            self.x.ge(&other.x)
        }
    }

    /// Trait to handle pairwise iteration over structs that can iter.
    pub trait Pairwise {
        // Where type is return of the function.
        fn pairwise_map<
            T: Send + Signed + Clone,
            F: Fn(Record, Record) -> Pair<Record, T> + Send + 'static + Sync,
        >(
            &self,
            function: F,
        ) -> Vec<Pair<Record, T>>;

        fn pairwise_map_option<
            T: Send + Signed + Clone,
            F: Fn(Record, Record) -> Option<Pair<Record, T>>
                + Send
                + 'static
                + Sync,
        >(
            &self,
            function: F,
        ) -> Vec<Pair<Record, T>>;

        fn pairwise_map_condensed_upper<
            T: Send + Signed + Clone,
            F: Fn(Record, Record) -> Pair<Record, T> + Send + 'static + Sync,
        >(
            &mut self,
            function: F,
        ) -> Vec<Pair<Record, T>>;
    }

    impl Pairwise for Vec<Record> {
        /// Map a function to all pairs in a paralel way

        fn pairwise_map<
            T: Send + Signed + Clone,
            F: Fn(Record, Record) -> Pair<Record, T> + Send + 'static + Sync,
        >(
            &self,
            function: F,
        ) -> Vec<Pair<Record, T>> {
            self.clone()
                .par_iter()
                .map(|x| {
                    self.par_iter().map(|y| function(x.clone(), y.clone()))
                })
                .flatten()
                .collect::<Vec<Pair<Record, T>>>()
        }

        /// Map a function to the pairwise combination of Fasta records in paralel.
        ///
        fn pairwise_map_option<
            T: Send + Signed + Clone,
            F: Fn(Record, Record) -> Option<Pair<Record, T>>
                + Send
                + 'static
                + Sync,
        >(
            &self,
            function: F,
        ) -> Vec<Pair<Record, T>> {
            let list = self
                .clone()
                .par_iter()
                .map(|x| {
                    self.par_iter()
                        .map(|y| function(x.clone(), y.clone()))
                        .filter_map(|x| x)
                })
                .flatten()
                .collect::<Vec<Pair<Record, T>>>();
            list
        }

        /// Map a function to the condensed upper triangle vector in paralel.
        ///
        /// First this calculates the combinations needed for the upper matrix and then
        /// applies the function in paralel to each of the combinations. Ignores the
        /// diagonal
        fn pairwise_map_condensed_upper<
            T: Send + Signed + Clone,
            F: Fn(Record, Record) -> Pair<Record, T> + Send + 'static + Sync,
        >(
            &mut self,
            function: F,
        ) -> Vec<Pair<Record, T>> {
            // first we get the combinations of the upper triangle.
            // then we can call a pairwise map on this.
            let mut to_calculate: Vec<(Record, Record)> = Vec::new();
            let mut skip = HashSet::new();
            // finding the upper triangle pairs.
            for x in self.clone() {
                let x = x.clone();
                for y in self.clone() {
                    if !skip.contains(&x.id().to_owned())
                        && !skip.contains(&y.id().to_owned())
                    {
                        // Remove the diagonal
                        if x != y {
                            to_calculate.push((x.clone(), y.clone()));
                        }
                    } else {
                        // No op
                    }
                }
                skip.insert(x.id().to_owned());
            }

            let upper_vector = to_calculate
                .par_iter()
                .map(|pair| function(pair.0.to_owned(), pair.1.to_owned()))
                .collect::<Vec<Pair<Record, T>>>();
            upper_vector
        }
    }
}

/// Remove the empty records from an array.
pub fn remove_empty(x: Vec<Record>) -> Vec<Record> {
    let x: Vec<Record> = x
        .par_iter()
        .filter(|x| !x.is_empty())
        .map(|x| x.to_owned())
        .collect();
    x
}

/// Module for the distance metrics
///
/// See the [mega](https://www.megasoftware.net/mega1_manual/Distance.html) manual for the
/// sources of these formulae.
///
/// Maybe for undefined distances I can use [std::f64::INFINITY].
///
/// See [emboss-tool](https://www.bioinformatics.nl/cgi-bin/emboss/distmat).
mod distances {
    use bio::alignment::distance;
    use bio::utils::TextSlice;

    use std::error::Error;

    /// The percentage of characters that is different between two strings.
    ///
    /// This is also known as the p-distance.
    pub fn percentage_difference(
        alpha: TextSlice,
        beta: TextSlice,
    ) -> Result<f64, Box<dyn Error>> {
        if alpha.len() != beta.len() {
            let string = format!(
                "Strings must be the same length: a: {} != b: {}",
                alpha.len(),
                beta.len()
            );

            return Err(string.into());
        }

        if alpha == beta {
            Ok(0.0)
        } else {
            let len = alpha.len();
            let distance = distance::hamming(alpha, beta) as f64;
            Ok(distance / len as f64)
        }
    }

    /// # Jukes Cantor distance
    ///
    /// Nucleotide substitution model. It assumes that all substitutions are equally likely. This
    /// is only defined for a [percentage_difference] of less than 0.75.
    ///
    /// # See
    /// * [mega](https://www.megasoftware.net/web_help_7/hc_jukes_cantor_distance.htm)
    /// * [rice uni](https://www.cs.rice.edu/~nakhleh/COMP571/Slides/Phylogenetics-DistanceMethods-Full.pdf)
    /// * [FUB](https://www.mi.fu-berlin.de/wiki/pub/ABI/PhylogeneticAnalsysis/probabilistic.pdf)
    pub fn jukes_cantor(
        alpha: TextSlice<'_>,
        beta: TextSlice<'_>,
    ) -> Result<f64, Box<dyn Error>> {
        // Calculate the proportion of differences (p)
        let p = percentage_difference(alpha, beta)?;

        // Check for edge case: p must be < 0.75
        if p >= 0.75 {
            let string = format!("Proportion of differences (p) is too high for Jukes-Cantor model: {p:.2} >= 0.75");

            return Err(string.into());
        }

        if p == 0.0 {
            Ok(0.0)
        } else {
            // Compute the distance using the Jukes-Cantor formula
            // d = - 3/4 ln (1 - 4/3 p)
            let x = 1.0 - (4.0 / 3.0 * p);
            let log = x.ln();
            let distance = -0.75 * log;
            Ok(distance)
        }
    }
}

/// JC distance between two fasta records `x`, and `y`.
///
/// This distance is only defined for two records of equal length.
pub fn fasta_distance_jukes_cantor(
    x: Record,
    y: Record,
) -> Option<Pair<Record, f64>> {
    let d = jukes_cantor(x.seq(), y.seq());
    match d {
        Ok(d) => Some(Pair::new(x, y, d)),
        Err(_) => None,
    }
}

/// JC distance between two fasta records `x`, and `y`.
///
/// This distance is only defined for two records of equal length and for a treshold of p > 75.
/// If this is not satisfied, [f64::InFINITY] is returned.
pub fn fasta_distance_jukes_cantor_number(
    x: Record,
    y: Record,
) -> Pair<Record, f64> {
    let d = jukes_cantor(x.seq(), y.seq());
    match d {
        Ok(d) => Pair::new(x, y, d),
        Err(e) => Pair::new(x, y, f64::INFINITY),
    }
}

/// Levenstein distance between two fasta records `x`, and `y`.
pub fn fasta_distance_levenshtein(
    x: Record,
    y: Record,
) -> Pair<Record, i64> {
    let d = levenshtein(x.seq(), y.seq()) as i64;
    Pair::new(x, y, d)
}

fn len_condensed_dis_mat(n: usize) -> usize {
    (n * (n - 1)) / 2
}

#[cfg(test)]
mod test_metrics {
    use approx_eq::assert_approx_eq;
    use crate::core::utils::distances;

    #[test]
    fn test_percentage_distance() {
        let a = b"ATGCTCGTAGCTGCTACGTC";
        let b = b"ATCCTCGAAGCAGCTACGAC";
        let d =
            distances::percentage_difference(a.as_ref(), b.as_ref()).unwrap();
        assert!(d == 0.20);

        let a = b"ATGCTAGC";
        let b = b"ATGCCAGT";
        let d =
            distances::percentage_difference(a.as_ref(), b.as_ref()).unwrap();
        assert!(d == 0.25);

        let a = b"ATGCTCGTAGCTGCTACGTCA";
        let b = b"ATCCTCGAAGCAGCTACGAC";
        let d = distances::percentage_difference(a.as_ref(), b.as_ref());
        assert!(d.is_err());
        if let Err(err) = d {
            assert_eq!(
                err.to_string(),
                "Strings must be the same length: a: 21 != b: 20"
            );
        }
    }

    /// This [page](http://www.insilicase.com/Web/JukesCantor.aspx) shows the correct results
    /// for the jukes cantor distance.
    #[test]
    fn test_jukes_cantor_distance() {
        // Correct calculations
        let a = b"ATGCTCGTAGCTGCTACGTC";
        let b = b"ATCCTCGAAGCAGCTACGAC";
        let d = distances::jukes_cantor(a.as_ref(), b.as_ref()).unwrap();
        assert_approx_eq!(d, 0.23261619622787);

        let a = b"ATGCTAGC";
        let b = b"ATGCCAGT";
        let d = distances::jukes_cantor(a.as_ref(), b.as_ref()).unwrap();
        assert_approx_eq!(d, 0.30409883108112);

        // Error condtions
        let a = b"ATGCTCGTAGCTGCTACGTCA";
        let b = b"ATCCTCGAAGCAGCTACGAC";
        let d = distances::jukes_cantor(a.as_ref(), b.as_ref());
        assert!(d.is_err());

        let a = b"TTTTAAC";
        let b = b"CCCCCCC";
        let d = distances::jukes_cantor(a.as_ref(), b.as_ref());
        assert!(d.is_err());
        if let Err(err) = d {
            assert_eq!(err.to_string(), "Proportion of differences (p) is too high for Jukes-Cantor model: 0.86 >= 0.75");
        }
    }
}

/// Recursive N-choose K
pub fn choose(
    n: f64,
    k: f64,
) -> f64 {
    if k == 0.0 {
        1.0
    } else {
        (n * choose(n - 1.0, k - 1.0)) / k
    }
}

#[cfg(test)]
mod test_traits {
    use approx_eq::assert_approx_eq;
    use bio::io::fasta::Record;

    use rayon::iter::IntoParallelRefIterator;
    use crate::core::io::read_fasta;
    use crate::core::utils::{choose, fasta_distance_jukes_cantor_number, remove_empty};
    use crate::core::utils::pairs::{Pair, Pairwise};

    #[test]
    fn test_choose() {
        assert_approx_eq!(choose(1.0, 2.0), 0.0);
        assert_approx_eq!(choose(10.0, 2.0), 45.0);
    }

    #[test]
    fn test_upper_distance() {
        use rayon::iter::ParallelIterator;

        fn test_n_matrix(
            file: &str,
            n: usize,
        ) {
            let file = file;
            match read_fasta(file) {
                Ok(x) => {
                    let len = x.len();
                    if len <= n {
                        panic!("too little test data {file} has {len} records: cannot take {n}")
                    }
                    let mut x: Vec<Record> =
                        remove_empty(x).iter().take(n).cloned().collect();
                    let expected_records = choose(n as f64, 2.0);
                    let distances = x.pairwise_map_condensed_upper(
                        fasta_distance_jukes_cantor_number,
                    );
                    let diagonal: Vec<Pair<Record, f64>> = distances
                        .par_iter()
                        .filter(|x| x.a == x.b)
                        .cloned()
                        .collect();
                    assert_eq!(diagonal.len(), 0, "There may be no diagonals");
                    assert_eq!(
                        expected_records,
                        distances.len() as f64,
                        "The matrix is false {distances:?}"
                    );
                }
                Err(e) => {
                    panic!("Error: {}", e)
                }
            };
        }

        test_n_matrix("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa", 2);
        test_n_matrix("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa", 3);
        test_n_matrix("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa", 4);
        test_n_matrix("resources/test/data/asv-listerria-taxon-Bacillales-Order.fasta.final_tree.fa", 15);
    }
}
