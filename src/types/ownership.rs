//! Ownership

/// Ownership
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Ownership {
    /// Owned by the current process
    Owned,
    /// Not owned on the current process. The two values are the process that owns it
    /// and its local index on that process
    Ghost(usize, usize),
}
