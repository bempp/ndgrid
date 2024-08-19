//! Point
use crate::{traits::Point as PointTrait, types::RealScalar};

/// A points
#[derive(Debug, Clone, Copy)]
pub struct Point<'a, T: RealScalar> {
    coordinates: &'a [T],
}

impl<'a, T: RealScalar> Point<'a, T> {
    /// Create new
    pub fn new(coordinates: &'a [T]) -> Self {
        Self { coordinates }
    }
}
impl<'a, T: RealScalar> PointTrait for Point<'a, T> {
    type T = T;

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
    points: Vec<&'a [T]>,
    index: usize,
}
impl<'a, T: RealScalar> PointIter<'a, T> {
    /// Create new
    pub fn new(points: Vec<&'a [T]>) -> Self {
        Self { points, index: 0 }
    }
}
impl<'a, T: RealScalar> Iterator for PointIter<'a, T> {
    type Item = Point<'a, T>;

    fn next(&mut self) -> Option<Point<'a, T>> {
        self.index += 1;
        if self.index <= self.points.len() {
            Some(Point::new(self.points[self.index - 1]))
        } else {
            None
        }
    }
}
