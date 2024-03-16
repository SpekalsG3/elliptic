use num_bigint::{BigInt, BigUint};
use num_traits::Zero;

// original from tutorial
#[deprecated(note="overflow when working close to boundary, use `u_xgcd`")]
pub fn xgcd (a: u128, b: u128) -> (i128, i128, u128) {
  let mut r = (a, b);
  let mut s = (1_i128, 0_i128);
  let mut t = (0_i128, 1_i128);

  while r.1 != 0 {
    let q = r.0 / r.1;
    r = (r.1, r.0 - q * r.1);
    s = (s.1, s.0 - q as i128 * s.1);
    t = (t.1, t.0 - q as i128 * t.1);
  }

  (s.0, t.0, r.0) // a, b, g
}

// version using unsigned arithmetics - no overflow
// source - https://jeffhurchalla.com/2018/10/13/implementing-the-extended-euclidean-algorithm-with-unsigned-inputs/
pub fn u_xgcd (a: BigUint, b: BigUint) -> (BigInt, BigInt, BigUint) {
  let mut xy1 = (BigInt::from(1_i128), BigInt::from(0_i128));

  let mut xy0 = (BigInt::from(0_i128), BigInt::from(1_i128));

  let mut a = (a, b);

  let mut q = BigUint::zero();

  while !a.1.is_zero() {
    let q_i = BigInt::from(q.clone());
    let x2 = xy0.0 - q_i.clone() * xy1.0.clone();
    let y2 = xy0.1 - q_i * xy1.1.clone();

    xy0 = xy1;
    xy1 = (x2, y2);

    // let a0 = a.0;
    // a.0 = a.1;
    // q = a0 / a.0;
    // a.1 = a0 - q * a.0;

    q = a.0.clone() / a.1.clone();
    a = (a.1.clone(), a.0 - q.clone() * a.1.clone())
  }

  (xy1.0, xy1.1, a.0.into())
}

#[cfg(test)]
mod tests {
  use crate::field::field::get_field_prime;
  use super::*;

  #[test]
  #[allow(deprecated)]
  fn test_xgcd () {
    assert_eq!(xgcd(10, 5), (0, 1, 5));
    assert_eq!(xgcd(240, 46), (-9, 47, 2));
    // assert_eq!(xgcd(3, 1 + BIG_PRIME), (1, 0, 1)); // fails
  }

  #[test]
  fn test_u_xgcd () {
    assert_eq!(u_xgcd(BigUint::from(10_u8), BigUint::from(5_u8)), (BigInt::from(0), BigInt::from(1), BigUint::from(5_u8)));
    assert_eq!(u_xgcd(BigUint::from(240_u8), BigUint::from(46_u8)), (BigInt::from(-9), BigInt::from(47), BigUint::from(2_u8)));
    assert_eq!(u_xgcd(BigUint::from(3_u8), get_field_prime()), (BigInt::from(90165965714076793378641578922350040406_u128), BigInt::from(-1), BigUint::from(1_u8)));
  }
}
