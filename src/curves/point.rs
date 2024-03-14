use std::ops::Neg;
use num_traits::One;
use crate::field::field::Field;
use crate::field::field_element::FieldElement;

#[derive(Debug, Clone, PartialEq)]
pub struct Point<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
    pub z: FieldElement<'a>,
}

impl<'a> Point<'a> {
    pub fn is_infinity(&self) -> bool {
        self.z.is_zero()
    }

    pub fn infinity(f: &'a Field) -> Self {
        Self {
            x: f.zero(),
            y: f.one(),
            z: f.zero(),
        }
    }
}

impl Neg for Point<'_> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: self.x,
            y: -self.y,
            z: self.z,
        }
    }
}
