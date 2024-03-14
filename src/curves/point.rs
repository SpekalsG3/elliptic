use num_bigint::BigUint;
use crate::field::field_element::FieldElement;

#[derive(Debug, Clone, PartialEq)]
pub struct Point<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
}

impl<'a> Point<'a> {
    pub fn is_infinity(&self) -> bool {
        self.y.is_zero()
    }
    pub fn is_identity(&self) -> bool {
        self.y.value == BigUint::from(1_u8)
    }
}

impl Neg for Point<'_> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: self.x,
            y: -self.y,
        }
    }
}
