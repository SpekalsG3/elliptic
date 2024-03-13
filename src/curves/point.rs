use crate::field::field_element::FieldElement;

#[derive(PartialEq)]
pub struct Point<'a> {
    pub x: FieldElement<'a>,
    pub y: FieldElement<'a>,
}

impl Point<'_> {
    pub fn is_infinity(&self) -> bool {
        self.y.is_zero()
    }
    pub fn is_identity(&self) -> bool {
        self.y.value == 1
    }
}
