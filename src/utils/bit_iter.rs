use std::usize;

#[derive(Debug, PartialOrd, PartialEq)]
pub struct BitIter<T>{
  done: bool,
  pub(crate) b: usize, // biggest
  pub(crate) l: usize, // least
  pub(crate) value: T,
}

impl BitIter<u32> {
  pub fn at_big(mut self) -> Self {
    let size = u32::BITS as usize;
    self.b = (0..size)
        .rposition(|i| (self.value >> i) & 1 == 1)
        .unwrap_or(self.l);
    self
  }
}
impl From<u32> for BitIter<u32> {
  fn from (value: u32) -> Self {
    BitIter {
      done: false,
      b: u32::BITS as usize - 1,
      l: 0,
      value,
    }
  }
}
impl Iterator for BitIter<u32> {
  type Item = bool;
  fn next (&mut self) -> Option<Self::Item> {
    if self.done {
      return None;
    }

    let bit = self.value >> self.b & 1;
    if self.b == self.l {
      self.done = true
    } else {
      self.b -= 1;
    }

    Some(bit == 1)
  }
}

impl DoubleEndedIterator for BitIter<u32> {
  fn next_back(&mut self) -> Option<Self::Item> {
    if self.done {
      return None;
    }

    let bit = self.value >> self.l & 1;
    if self.b == self.l {
      self.done = true
    } else {
      self.l += 1;
    }

    Some(bit == 1)
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn all () {
    let mut s = BitIter::from(11);

    for _ in 0..(u32::BITS - 4) {
      assert_eq!(s.next(), Some(false));
    }
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), Some(false));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), None);
  }

  #[test]
  fn next () {
    let mut s = BitIter::from(11).at_big();

    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), Some(false));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), Some(true));
    assert_eq!(s.next(), None);
  }

  #[test]
  fn back () {
    let mut s = BitIter::from(11).at_big();

    assert_eq!(s.next_back(), Some(true));
    assert_eq!(s.next_back(), Some(true));
    assert_eq!(s.next_back(), Some(false));
    assert_eq!(s.next_back(), Some(true));
    assert_eq!(s.next_back(), None);
  }

  #[test]
  fn both () {
    let mut s = BitIter::from(11).at_big();

    assert_eq!(s.next(),      Some(true));
    assert_eq!(s.next_back(), Some(true));
    assert_eq!(s.next_back(), Some(true));
    assert_eq!(s.next(),      Some(false));
    assert_eq!(s.next_back(), None);
    assert_eq!(s.next(),      None);
  }

  #[test]
  fn zero () {
    let mut s = BitIter::from(0).at_big(); // 0

    assert_eq!(s.next(), Some(false));
    assert_eq!(s.next(), None);
  }
}
