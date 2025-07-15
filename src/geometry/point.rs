//! Point
use crate::{traits::Point as PointTrait, types::RealScalar};

/// A points
#[derive(Debug, Clone, Copy)]
pub struct Point<'a, T: RealScalar> {
    index: usize,
    coordinates: &'a [T],
}

impl<'a, T: RealScalar> Point<'a, T> {
    /// Create new
    pub fn new(index: usize, coordinates: &'a [T]) -> Self {
        Self { index, coordinates }
    }
}
impl<T: RealScalar> PointTrait for Point<'_, T> {
    type T = T;

    fn index(&self) -> usize {
        self.index
    }

    fn dim(&self) -> usize {
        self.coordinates.len()
    }

    fn coords(&self, data: &mut [T]) {
        data.copy_from_slice(self.coordinates);
    }
}

/// Iterator over points
#[derive(Debug)]
pub struct PointIter<'a, T: RealScalar> {
    points: Vec<(usize, &'a [T])>,
    index: usize,
}
impl<'a, T: RealScalar> PointIter<'a, T> {
    /// Create new
    pub fn new(points: Vec<(usize, &'a [T])>) -> Self {
        Self { points, index: 0 }
    }
}
impl<'a, T: RealScalar> Iterator for PointIter<'a, T> {
    type Item = Point<'a, T>;

    fn next(&mut self) -> Option<Point<'a, T>> {
        self.index += 1;
        if self.index <= self.points.len() {
            Some(Point::new(
                self.points[self.index - 1].0,
                self.points[self.index - 1].1,
            ))
        } else {
            None
        }
    }
}
