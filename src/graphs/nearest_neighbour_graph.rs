// use anyhow::Result;
// use ndarray::Array2;
// use hnsw_rs::{hnsw, dist};

// /// constructs a nearest neighbour graph from the input table
// /// and return the row indices of the contigs that are determined
// /// to be too distant from their nearest neighbours
// pub fn find_distant_elements(m: &clap::ArgMatches, input_array: Array2<f32>) -> Result<Vec<usize>> {
//     let n_neighbours = *m.get_one::<usize>("n-neighbours").unwrap();
//     let mut distant_elements = Vec::new();
    
//     let mut graph = hnsw::Hnsw::new(
//         n_neighbours, 
//         16, 
//         200, 
//         hnsw_rs::DistanceType::Dot,

//     )?;
//     graph.build(16, 200, hnsw_rs::DistanceType::Dot)?;
//     for (i, row) in input_array.rows().into_iter().enumerate() {
//         let nearest_neighbour = graph.search_knn(&row, 1, hnsw_rs::DistanceType::Dot)?;
//         if nearest_neighbour[0].1 > 0.5 {
//             distant_elements.push(i);
//         }
//     }

//     Ok(distant_elements)
// } 