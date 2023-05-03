use anyhow::Result;
use ndarray::{Array2, Data};
use linfa_nn::{CommonNearestNeighbour, distance::Distance, NearestNeighbour};
use linfa::Float;

pub fn nearest_neighbour_filter(coverage_array: &Array2<f32>, kmer_array: &Array2<f32>) -> Vec<usize> {
    Vec::new()
}

fn kd_tree_filter<F: Float, D: Distance<F> + Clone>(array: &Array2<F>, function: D, threshold: F) -> Result<Vec<usize>> {
    let mut disconnected_points = Vec::new();
    let nn = CommonNearestNeighbour::KdTree.from_batch(array, function.clone())?;
    for (i, row) in array.rows().into_iter().enumerate() {
        let nearest_neighbour = nn.k_nearest(row.view(), 1)?;
        // calculate the distance between the point and its nearest neighbour
        let distance = function.distance(nearest_neighbour[0].0.view(), row.view());
        if distance > threshold {
            disconnected_points.push(i);
        }
    }
    Ok(disconnected_points)
}

const nearest_neighbors: &str = r#"
import numpy as np
import warnings
from sklearn.manifold import SpectralEmbedding
from pynndescent import NNDescent

from umap.umap_ import nearest_neighbors
from umap.umap_ import fuzzy_simplicial_set
from umap.umap_ import find_ab_params

mimport warnings
warnings.filterwarnings('ignore')

#Bootstrap function to do hyperparam search
def full_umap(X, metric, min_dist, dim, init_nn, to_add_nn, connectivity, rand_state, y):
  
  random_state = check_random_state(rand_state)
  
  a, b = find_ab_params(1.0, min_dist)
  
  knn_indices, knn_dists, knn_search_index = nearest_neighbors(
    X,
    n_neighbors=init_nn,
    metric = metric,
    metric_kwds = {},
    angular=False,
    random_state = random_state,
    low_memory=True,
    use_pynndescent=True,
    n_jobs=1,
    verbose=False,
  )
  
  
  if connectivity == "org":
    new_knn_indices, new_knn_dists = knn_indices, knn_dists
  else :
    connected_mnn, mutual_nn, new_knn_indices, new_knn_dists  = mutual_nn_nearest(knn_indices, knn_dists, init_nn, to_add_nn, connectivity=connectivity)
  
  
  
  P, sigmas, rhos = fuzzy_simplicial_set(
    X = X,
    n_neighbors = to_add_nn,
    metric = metric,
    random_state = random_state,
    knn_indices= new_knn_indices,
    knn_dists = new_knn_dists,
  )

  embeddings, aux_data = simplicial_set_embedding(
    data = X,
    graph = P,
    n_components = dim,
    initial_alpha = 1,
    a = a,
    b = b,
    gamma = 1.0,
    negative_sample_rate = 5,
    n_epochs = 200,
    init = "spectral",
    random_state = check_random_state(rand_state),
    metric = metric,
    metric_kwds = {},
    densmap = False,
    densmap_kwds = {},
    output_dens = False,
    output_metric= dist.named_distances_with_gradients["euclidean"],
    output_metric_kwds={},
    euclidean_output=True,
    parallel=False,
    verbose=False,
  )
  embeddings = np.nan_to_num(embeddings)
  return embeddings
"#;

const minimum_spanning_tree: &str = r#"

#Calculate min spanning tree

def min_spanning_tree(knn_indices, knn_dists, n_neighbors, threshold):
  
  rows = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
  cols = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.int32)
  vals = np.zeros(knn_indices.shape[0] * n_neighbors, dtype=np.float32)
  
  pos = 0
  for i, indices in enumerate(knn_indices):
    for j, index in enumerate(indices[:threshold]):
      if index == -1:
        continue
      rows[pos] = i 
      cols[pos] = index
      vals[pos] = knn_dists[i][j]
      pos += 1
  
  matrix = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(knn_indices.shape[0], knn_indices.shape[0]))
  Tcsr = scipy.sparse.csgraph.minimum_spanning_tree(matrix)
  
  Tcsr = scipy.sparse.coo_matrix(Tcsr)
  weights_tuples = zip(Tcsr.row, Tcsr.col, Tcsr.data)
  

  sorted_weights_tuples = sorted(weights_tuples, key=lambda tup: tup[2])
  
  return sorted_weights_tuples 

"#;