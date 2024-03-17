use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use num_bigint::{BigInt, BigUint};
use num_traits::Zero;
use crate::field::field_element::FieldElement;
use crate::utils::xgcd::u_xgcd;

#[cfg(test)]
pub fn get_field_prime() -> BigUint {
  // 270497897142230380135924736767050121217
  BigUint::from(1_u128 + 407 * (1 << 119))
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Field {
  pub order: BigUint,
}

impl From<&str> for Field {
  fn from(value: &str) -> Self {
    Self {
      order: value.parse().unwrap(),
    }
  }
}

impl Display for Field {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.order)
  }
}

impl<'a> Field {
  pub fn new (order: BigUint) -> Field {
    // assert_eq!(order, FIELD_PRIME, "Only 1+407*2^119 currently implemented");
    Field {
      order,
    }
  }

  pub fn get(&self, v: BigUint) -> FieldElement {
    FieldElement::new(self, v % self.order.clone())
  }

  pub fn zero (&'a self) -> FieldElement<'a> {
    FieldElement {
      field: &self,
      value: BigUint::zero(),
    }
  }

  pub fn one (&'a self) -> FieldElement<'a> {
    FieldElement {
      field: &self,
      value: BigUint::from(1_u8),
    }
  }

  pub fn sample (&self, bytes: &[u8]) -> FieldElement<'_> {
    let res = bytes
        .iter()
        .fold(BigUint::zero(), |acc, b| {
          (acc << 8) ^ BigUint::from(*b)
        });

    FieldElement::new(self, res % self.order.clone())
  }

  pub(crate) fn sub_mod (&self, a: BigUint, b: BigUint) -> BigUint {
    match a.cmp(&b) {
      Ordering::Greater => a - b,
      Ordering::Equal => BigUint::zero(),
      Ordering::Less => self.neg_mod(b - a),
    }
  }

  pub(crate) fn add_mod (&self, a: BigUint, b: BigUint) -> BigUint {
    self.sub_mod(a, self.order.clone() - b)
  }

  pub(crate) fn mul_mod (&self, a: BigUint, b: BigUint) -> BigUint {
    let big_zero = BigUint::zero();
    let big_one = BigUint::from(1_u8);

    let mut res = BigUint::zero();

    let mut a = a;
    let mut b = b;

    while b > big_zero {
      if b.clone() & big_one.clone() == big_one {
        res = self.add_mod(res, a.clone());
      }

      a = self.add_mod(a.clone(), a);

      b /= BigUint::from(2_u8);
    }

    // // source - https://www.youtube.com/watch?v=9hSmQtL49g4&ab_channel=DG
    // // doesnt work
    // while a != 0 {
    //   if a & 1 == 1  {
    //     res ^= b;
    //   }
    //
    //   a >>= 1;
    //   b <<= 1;
    //
    //   if degree(b) == self.degree {
    //     b ^= self.prime;
    //   }
    // }

    res
  }

  pub(crate) fn neg_mod (&self, a: BigUint) -> BigUint {
    if a.is_zero() {
      a
    } else {
      self.order.clone() - a
    }
  }

  // inverse of `x` is `x ** -1 = 1/x` so that `x` multiplied by inversed `x` is `1`
  pub(crate) fn inv (&self, a: BigUint) -> BigUint {
    let (a, _, _) = u_xgcd(a, self.order.clone());

    // because a can be negative
    match a.cmp(&BigInt::zero()) {
      Ordering::Greater => a.to_biguint().unwrap(),
      Ordering::Equal => BigUint::zero(),
      Ordering::Less => self.neg_mod((-a).to_biguint().unwrap()),
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn mul () {
    let field = Field::new(get_field_prime());

    assert_eq!(field.mul_mod(BigUint::from(2_u128), BigUint::from(3_u128)), BigUint::from(6_u8));
    assert_eq!(field.mul_mod(get_field_prime(), BigUint::from(3_u128)), BigUint::from(0_u8));
    assert_eq!(field.mul_mod(get_field_prime() - BigUint::from(1_u8), BigUint::from(3_u128)), get_field_prime() - BigUint::from(3_u8));
  }

  #[test]
  fn neg () {
    let field = Field::new(get_field_prime());
    assert_eq!(field.neg_mod(BigUint::from(256_u128)), BigUint::from(270497897142230380135924736767050120961_u128));
    let field = Field::new(BigUint::from(100_u128));
    assert_eq!(
      field.add_mod(BigUint::from(20_u128), field.neg_mod(BigUint::from(20_u128))),
      BigUint::from(0_u8),
    );
    assert_eq!(
      field.add_mod(BigUint::from(20_u128), field.neg_mod(BigUint::from(19_u128))),
      BigUint::from(1_u8),
    );
  }
}
