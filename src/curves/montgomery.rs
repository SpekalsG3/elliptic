use std::str::FromStr;
use num_bigint::BigUint;
use crate::curves::point::Point;
use crate::field::field_element::FieldElement;
use crate::utils::bit_iter::BitIter;
use crate::utils::s_tonelli::tonelli_shanks;

pub const MONTGOMERY_A: u128 = 486662;
pub const MONTGOMERY_B: u128 = 1;
// fn get_montgomery_order() -> BigUint {
//     // (1 << 252) + 27742317777372353535851937790883648493
//     BigUint::from_str("7237005577332262213973186563042994240857116359379907606001950938285454250989").unwrap()
// }
fn get_montgomery_prime() -> BigUint {
    // (1_u128 << 255_u128) - 19
    BigUint::from_str("57896044618658097711785492504343953926634992332820282019728792003956564819949").unwrap()
}

// B * y^2 = x^3 + A * x^2 + x
pub struct MontgomeryCurve<'a> {
    a: FieldElement<'a>,
    b: FieldElement<'a>,
}

// https://www.hyperelliptic.org/EFD/g1p/auto-montgom.html
impl<'a> MontgomeryCurve<'a> {
    pub fn new(
        a: FieldElement<'a>,
        b: FieldElement<'a>,
    ) -> Self {
        assert_eq!(a.field.order, b.field.order, "has to be in the same field");

        Self {
            a,
            b,
        }
    }

    pub fn generator(
        &self,
    ) -> Point<'a> {
        assert_eq!(self.a.value, MONTGOMERY_A.into(), "unknown curve");
        assert_eq!(self.b.value, MONTGOMERY_B.into(), "unknown curve");
        assert_eq!(self.a.field.order, get_montgomery_prime(), "unknown field");

        let x = FieldElement::new(self.a.field, BigUint::from(9_u8));
        Point {
            x: x.clone(),
            y: tonelli_shanks((x.clone() ^ BigUint::from(3_u8)) + self.a.clone() * (x.clone() ^ BigUint::from(2_u8)) + x).expect("should be solvable").0,
        }
    }

    pub fn point_double(
        &self,
        p1: &Point<'a>,
    ) -> Point<'a> {
        let one = self.a.field.one();
        let two = FieldElement::new(self.a.field, BigUint::from(2_u8));
        let three = FieldElement::new(self.a.field, BigUint::from(3_u8));

        let slope = (three * (p1.x.clone() ^ BigUint::from(2_u8)) + two.clone() * self.a.clone() * p1.x.clone() + one) / (two.clone() * self.b.clone() * p1.y.clone());
        let x = self.b.clone() * (slope.clone() ^ BigUint::from(2_u8)) - self.a.clone() - p1.x.clone() - p1.x.clone();
        let y = (two * p1.x.clone() + p1.x.clone() + self.a.clone()) * slope.clone() - self.b.clone() * (slope ^ BigUint::from(3_u8)) - p1.y.clone();
        Point { x, y }
    }

    pub fn point_add(
        &self,
        p1: &Point<'a>,
        p2: &Point<'a>,
    ) -> Point<'a> {
        if p1.x == p2.x {
            if p1.y == -p2.y.clone() {
                return Point {
                    x: p1.x.clone(),
                    y: p1.x.field.zero(),
                }
            }
            return self.point_double(p1);
        }

        let two = FieldElement::new(self.a.field, BigUint::from(2_u8));

        let slope = (p2.y.clone() - p1.y.clone()) / (p2.x.clone() - p1.x.clone());
        let x = self.b.clone() * (slope.clone() ^ BigUint::from(2_u8)) - self.a.clone() - p1.x.clone() - p2.x.clone();
        let y = (two * p1.x.clone() + p2.x.clone() + self.a.clone()) * slope.clone() - self.b.clone() * (slope ^ BigUint::from(3_u8)) - p1.y.clone();
        Point { x, y }
    }

    pub fn ladder(
        &self,
        n: FieldElement<'a>,
        p1: Point<'a>,
    ) -> Point<'a> {
        let mut r0 = Point { x: p1.x.field.zero(), y: p1.x.field.zero() };
        let mut r1 = p1;
        for digit in n.value.iter_u32_digits() {
            for b in BitIter::new(u128::BITS as usize, digit) {
                if b {
                    r0 = self.point_add(&r0, &r1);
                    r1 = self.point_double(&r1);
                } else {
                    r1 = self.point_add(&r0, &r1);
                    r0 = self.point_double(&r0);
                }
            }
        }
        r0
    }
}
