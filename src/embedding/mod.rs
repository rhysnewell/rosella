pub mod features;
pub mod intersect;
pub mod knn;
pub mod layout;
pub mod manifold;
pub mod metrics;
pub mod quality;
pub mod spectral;
pub mod umap;

pub type Graph = sprs::CsMatI<f32, u32, usize>;
