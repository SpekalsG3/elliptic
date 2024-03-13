use crate::curves::point::Point;
use crate::field::field_element::FieldElement;
use crate::utils::bit_iter::BitIter;
use crate::utils::s_tonelli::tonelli_shanks;

pub const MONTGOMERY_A: u128 = 486662;
pub const MONTGOMERY_B: u128 = 1;
// pub const MONTGOMERY_ORDER: u128 = (1 << 252) + 27742317777372353535851937790883648493;
// pub const MONTGOMERY_PRIME: u128 = (1_u128 << 255_u128) - 19;
pub const MONTGOMERY_PRIME: u128 = 3;

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
        assert_eq!(self.a.value, MONTGOMERY_A, "unknown curve");
        assert_eq!(self.b.value, MONTGOMERY_B, "unknown curve");
        assert_eq!(self.a.field.order, MONTGOMERY_PRIME, "unknown field");

        let x = FieldElement::new(self.a.field, 9);
        Point {
            x,
            y: tonelli_shanks((x ^ 3_usize) + self.a * (x ^ 2_usize) + x).expect("should be solvable").0,
        }
    }

    pub fn point_double(
        &self,
        p1: &Point<'a>,
    ) -> Point<'a> {
        let one = self.a.field.one();
        let two = FieldElement::new(self.a.field, 2);
        let three = FieldElement::new(self.a.field, 3);

        let slope = (three * (p1.x ^ 2_usize) + two * self.a * p1.x + one) / (two * self.b * p1.y);
        let x = self.b * (slope ^ 2_usize) - self.a - p1.x - p1.x;
        let y = (two * p1.x + p1.x + self.a) * slope - self.b * (slope ^ 3_usize) - p1.y;
        Point { x, y }
    }

    pub fn point_add(
        &self,
        p1: &Point<'a>,
        p2: &Point<'a>,
    ) -> Point<'a> {
        if p1.x == p2.x {
            if p1.y == -p2.y {
                return Point {
                    x: p1.x,
                    y: p1.x.field.zero(),
                }
            }
            return self.point_double(p1);
        }

        let two = FieldElement::new(self.a.field, 2);

        let slope = (p2.y - p1.y) / (p2.x - p1.x);
        let x = self.b * (slope ^ 2_usize) - self.a - p1.x - p2.x;
        let y = (two * p1.x + p2.x + self.a) * slope - self.b * (slope ^ 3_usize) - p1.y;
        Point { x, y }
    }

    pub fn ladder(
        &self,
        n: FieldElement<'a>,
        p1: Point<'a>,
    ) -> Point<'a> {
        let mut r0 = Point { x: p1.x.field.zero(), y: p1.x.field.zero() };
        let mut r1 = p1;
        for b in BitIter::new(u128::BITS as usize, n.value) {
            if b {
                r0 = self.point_add(&r0, &r1);
                r1 = self.point_double(&r1);
            } else {
                r1 = self.point_add(&r0, &r1);
                r0 = self.point_double(&r0);
            }
        }
        r0
    }
}
