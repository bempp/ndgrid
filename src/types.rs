//! Types

mod ownership;
pub use ownership::Ownership;
mod connectivity;
pub use connectivity::CellLocalIndexPair;

use num::Float;
use rlst::{LinAlg, RlstScalar};

/// A real scalar
pub trait RealScalar: Float + LinAlg + RlstScalar<Real = Self> {}

/// A simple integer matrix type for storing indices
#[derive(Debug)]
pub struct IntegerArray2 {
    data: Vec<usize>,
    dim: [usize; 2],
}

impl IntegerArray2 {
    /// Create a new integer array
    pub fn new(dim: [usize; 2]) -> Self {
        let nelems = dim.iter().product();
        Self {
            data: vec![0; nelems],
            dim,
        }
    }

    /// Create a new integer array from a slice
    pub fn new_from_slice(data: &[usize], dim: [usize; 2]) -> Self {
        let nelems = dim.iter().product();
        assert_eq!(data.len(), nelems);
        Self {
            data: data.to_vec(),
            dim,
        }
    }

    /// The shape of the array
    pub fn dim(&self) -> [usize; 2] {
        self.dim
    }

    /// Column iterator
    pub fn col_iter(&self) -> ColIter<'_> {
        ColIter {
            arr: self,
            index: 0,
        }
    }

    /// Mutable column iterator
    pub fn col_iter_mut(&mut self) -> ColIterMut<'_> {
        ColIterMut {
            arr: self,
            index: 0,
        }
    }
}

impl std::ops::Index<[usize; 2]> for IntegerArray2 {
    type Output = usize;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.data[self.dim[0] * index[1] + index[0]]
    }
}

impl std::ops::IndexMut<[usize; 2]> for IntegerArray2 {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.data[self.dim[0] * index[1] + index[0]]
    }
}

/// Column iterator
pub struct ColIter<'a> {
    arr: &'a IntegerArray2,
    index: usize,
}

impl<'a> std::iter::Iterator for ColIter<'a> {
    type Item = &'a [usize];

    fn next(&mut self) -> Option<Self::Item> {
        let nrows = self.arr.dim()[0];
        let index = self.index;
        self.index += 1;
        if index < self.arr.dim[1] {
            Some(&self.arr.data[index * nrows..(1 + index) * nrows])
        } else {
            None
        }
    }
}

/// Mutable column iterator
pub struct ColIterMut<'a> {
    arr: &'a mut IntegerArray2,
    index: usize,
}

impl<'a> std::iter::Iterator for ColIterMut<'a> {
    type Item = &'a mut [usize];

    fn next(&mut self) -> Option<Self::Item> {
        let nrows = self.arr.dim()[0];
        let index = self.index;
        self.index += 1;
        if index < self.arr.dim[1] {
            Some(unsafe {
                std::mem::transmute(&mut self.arr.data[index * nrows..(1 + index) * nrows])
            })
        } else {
            None
        }
    }
}
