use num_bigint::BigUint;
use crate::curves::point::Point;
use crate::field::field_element::FieldElement;

// y^2 = x^3 + Ax + B
pub struct WeierstrassCurve<'a> {
    a: FieldElement<'a>,
    #[allow(unused)]
    b: FieldElement<'a>,
}

// https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
impl<'a> WeierstrassCurve<'a> {
    pub fn new (
        a: FieldElement<'a>,
        b: FieldElement<'a>,
    ) -> Self {
        assert_eq!(a.field, b.field, "should be in the same field");
        Self {
            a,
            b,
        }
    }

    pub fn point_double(
        &self,
        p1: Point<'a>,
    ) -> Point<'a> {
        let three = FieldElement::new(p1.x.field, BigUint::from(3_u8));
        let slope = (three * (p1.x.clone() ^ BigUint::from(2_usize)) + self.a.clone()) / (p1.y.clone() + p1.y.clone());
        let x = (slope.clone() ^ BigUint::from(2_usize)) - (p1.x.clone() + p1.x.clone());
        Point {
            x: x.clone(),
            y: slope * (p1.x - x) - p1.y,
        }
    }

    pub fn point_add(
        &self,
        p1: Point<'a>,
        p2: Point<'a>,
    ) -> Point<'a> {
        if p1.x == p2.x {
            if p1.y == -p2.y {
                return Point {
                    x: p1.x.clone(),
                    y: p1.x.field.zero(),
                }
            }
            return self.point_double(p1);
        }
        let slope = (p2.y - p1.y) / (p2.x.clone() - p1.x.clone());
        let x = (slope.clone() ^ BigUint::from(2_usize)) - p1.x.clone() - p2.x;
        Point {
            x: x.clone(),
            y: slope * (p1.x - x),
        }
    }
}
