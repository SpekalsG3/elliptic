use std::fmt::{Debug, Formatter};
use std::ops::{Add, BitXor, Div, Mul, Neg, Sub};
use num_bigint::BigUint;
use num_traits::Zero;
use crate::field::field::Field;

#[derive(Clone, PartialEq, PartialOrd)]
pub struct FieldElement<'a> {
  pub field: &'a Field,
  pub value: BigUint,
}

impl Debug for FieldElement<'_> {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.value)
  }
}

impl<'a> FieldElement<'a> {
  pub fn unstringify(str: &str, field: &'a Field) -> Self {
    Self {
      field,
      value: str.parse().unwrap(),
    }
  }

  pub fn new(field: &'a Field, value: BigUint) -> Self {
    Self {
      field,
      value,
    }
  }

  pub fn inverse (&self) -> FieldElement<'a> {
    FieldElement {
      field: self.field,
      value: self.field.inv(self.value.clone()),
    }
  }

  pub fn is_zero (&self) -> bool {
    self.value.is_zero()
  }
}

impl<'a> Add for FieldElement<'a> {
  type Output = Self;
  fn add (self, rhs: Self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.add_mod(self.value, rhs.value),
    }
  }
}

impl<'a> Sub for FieldElement<'a> {
  type Output = Self;
  fn sub (self, rhs: Self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.sub_mod(self.value, rhs.value),
    }
  }
}

impl<'a> Mul for FieldElement<'a> {
  type Output = Self;
  fn mul (self, rhs: Self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.mul_mod(self.value, rhs.value),
    }
  }
}

impl<'a> Div for FieldElement<'a> {
  type Output = FieldElement<'a>;
  fn div (self, rhs: Self) -> Self::Output {
    assert_ne!(rhs.value, BigUint::zero(), "divide by zero");

    self * rhs.inverse()
  }
}

impl<'a> Neg for FieldElement<'a> {
  type Output = FieldElement<'a>;
  fn neg (self) -> Self::Output {
    FieldElement {
      field: self.field,
      value: self.field.neg_mod(self.value),
    }
  }
}

impl<'a> ToString for FieldElement<'a> {
  fn to_string (&self) -> String {
    self.value.to_string()
  }
}

// todo tutorial uses `BitXor` for power, replaces later
impl<'a> BitXor<BigUint> for FieldElement<'a> {
  type Output = Self;

  fn bitxor (self, exponent: BigUint) -> Self::Output {
    Self::new(self.field, self.value.modpow(&exponent, &self.field.order))
  }
}

#[cfg(test)]
mod tests {
  use crate::field::field::get_field_prime;
  use super::*;

  #[test]
  fn mul () {
    let field = Field::new(get_field_prime());
    assert_eq!(
      FieldElement::new(&field, BigUint::from(49789714223038013592473676705012096123_u128)) * FieldElement::new(&field, BigUint::from(6534789852937546098347957826345234_u128)),
      FieldElement::new(&field, BigUint::from(105250150227149389100670877502232671566_u128)),
    );

    let el_1 = FieldElement::new(&field, BigUint::from(8_u128));
    let el_3 = FieldElement::new(&field, BigUint::from(12_u128));
    assert_eq!(el_1 * el_3, FieldElement::new(&field, BigUint::from(96_u128)));

    let el_1 = FieldElement::new(&field, BigUint::from(3_u128));
    let el_3 = FieldElement::new(&field, BigUint::from(270497897142230380135924736767050121215_u128));
    assert_eq!(el_1 * el_3, FieldElement::new(&field, BigUint::from(270497897142230380135924736767050121211_u128)));
  }

  #[test]
  fn div () {
    let field = Field::new(get_field_prime());
    // assert_eq!(
    //   FieldElement::new(&field, BigUint::from(74658620945386735627456854792784352353_u128)) / FieldElement::new(&field, BigUint::from(85408008396924667383611388730472331217_u128)),
    //   FieldElement::new(&field, BigUint::from(120557879365253444230411244907275635216_u128))
    // );

    let el_1 = FieldElement::new(&field, BigUint::from(12_u128));
    let el_2 = FieldElement::new(&field, BigUint::from(4_u128));
    assert_eq!(el_1 / el_2, FieldElement::new(&field, BigUint::from(3_u128)));

    let el_1 = FieldElement::new(&field, BigUint::from(270497897142230380135924736767050121215_u128));
    let el_2 = FieldElement::new(&field, BigUint::from(5_u128));
    assert_eq!(el_1 / el_2, FieldElement::new(&field, BigUint::from(54099579428446076027184947353410024243_u128)));

    let el_1 = FieldElement::new(&field, BigUint::from(270497897142230380135924736767050121215_u128));
    let el_2 = FieldElement::new(&field, BigUint::from(5_u128));
    assert_eq!(el_1 / el_2, FieldElement::new(&field, BigUint::from(54099579428446076027184947353410024243_u128)));

    assert_eq!(
      FieldElement::new(&field, BigUint::from(5012096123_u128)) / FieldElement::new(&field, BigUint::from(6534789852937546098347957826345234_u128)),
      FieldElement::new(&field, BigUint::from(109071144973379706934869779239844248849_u128)),
    );

    let field = Field::new(BigUint::from(8_u128));
    let el_1 = FieldElement::new(&field, BigUint::from(2_u128));
    let el_2 = FieldElement::new(&field, BigUint::from(7_u128));
    assert_eq!(el_1 / el_2, FieldElement::new(
      &field,
      BigUint::from(6_u8), // because 6 * 7 = 2 (mod 8)
    ));
  }

  #[test]
  fn add () {
    let field = Field::new(get_field_prime());
    assert_eq!(
      FieldElement::new(&field, BigUint::from(270497897142230380135924736767050120961_u128)) + FieldElement::new(&field, BigUint::from(300_u128)),
      FieldElement::new(&field, BigUint::from(44_u128)),
    );

    let field = Field::new(BigUint::from(100_u128));
    assert_eq!(
      FieldElement::new(&field, BigUint::from(20_u128)) + FieldElement::new(&field, BigUint::from(20_u128)),
      FieldElement::new(&field, BigUint::from(40_u128)),
    );
    assert_eq!(
      FieldElement::new(&field, BigUint::from(20_u128)) + FieldElement::new(&field, BigUint::from(19_u128)).neg(),
      field.one(),
    );
    assert_eq!(
      FieldElement::new(&field, BigUint::from(80_u128)) + FieldElement::new(&field, BigUint::from(21_u128)),
      field.one(),
    );
  }

  #[test]
  fn sub () {
    let field = Field::new(get_field_prime());
    assert_eq!(
      FieldElement::new(&field, BigUint::from(44_u128)) - FieldElement::new(&field, BigUint::from(200_u128)),
      FieldElement::new(&field, BigUint::from(270497897142230380135924736767050121061_u128)),
    );

    let field = Field::new(BigUint::from(100_u128));
    assert_eq!(
      FieldElement::new(&field, BigUint::from(20_u128)) - FieldElement::new(&field, BigUint::from(20_u128)),
      field.zero(),
    );
    assert_eq!(
      FieldElement::new(&field, BigUint::from(20_u128)) - FieldElement::new(&field, BigUint::from(19_u128)),
      field.one(),
    );
    assert_eq!(
      FieldElement::new(&field, BigUint::from(20_u128)) - FieldElement::new(&field, BigUint::from(21_u128)),
      field.one().neg(),
    );
  }

  #[test]
  fn neg () {
    let field = Field::new(get_field_prime());
    assert_eq!(
      FieldElement::new(&field, BigUint::from(6534789852937546098_u128)).neg(),
      FieldElement::new(&field, BigUint::from(270497897142230380129389946914112575119_u128)),
    );

    let field = Field::new(BigUint::from(100_u128));
    assert_eq!(
      FieldElement::new(&field, BigUint::from(1_u128)).neg(),
      FieldElement::new(&field, BigUint::from(99_u128)),
    );
    assert_eq!(
      FieldElement::new(&field, BigUint::from(20_u128)).neg(),
      FieldElement::new(&field, BigUint::from(80_u128)),
    );
  }

  #[test]
  fn pow () {
    let field = Field::new(get_field_prime());

    assert_eq!(FieldElement::new(&field, BigUint::from(6534789852937546098_u128)) ^ BigUint::from(501209126122_usize), FieldElement::new(&field, BigUint::from(256557788041265930815463337858691703671_u128)));

    let el = FieldElement::new(&field, BigUint::from(15_u8));
    assert_eq!(el ^ BigUint::from(4_usize), FieldElement::new(&field, BigUint::from(50625_u128)));

    let el = FieldElement::new(&field, BigUint::from(270497897142230380135_u128));
    assert_eq!(el ^ BigUint::from(8_usize), FieldElement::new(&field, BigUint::from(79016866124691016201920330826259043252_u128)));
  }

  // fn bitxor () {
  //   let field = Field::new(FIELD_PRIME);
  //
  //   let el = FieldElement {
  //     field: &field,
  //     value: 0b110011010, // 410
  //   };
  //   assert_eq!(el ^ 0b101101, FieldElement { // 45
  //     field: &field,
  //     value: 0b110110111, // 439
  //   });
  // }
}
