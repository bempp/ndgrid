//! Ownership

/// Ownership
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Ownership {
    /// Owned by the current process
    Owned,
    /// Ghost on the current process. The two values are the process that owns this and the local index on that process
    Ghost(usize, usize),
}
