use crate::curves::point::Point;
use crate::field::field_element::FieldElement;

// y^2 = x^3 + Ax + B
pub struct WeierstrassCurve<'a> {
    a: FieldElement<'a>,
    b: FieldElement<'a>,
}

// https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
impl<'a> WeierstrassCurve<'a> {
    pub fn new (
        a: FieldElement<'a>,
        b: FieldElement<'a>,
    ) -> Self {
        assert_eq!(a.field, b.field, "should be in the same field");
        Self { a, b }
    }

    pub fn evaluate (
        &self,
        x: FieldElement<'a>,
    ) -> FieldElement<'a> {
        x * x * x + self.a * x + self.b
    }

    pub fn point_double(
        &self,
        p1: Point<'a>,
    ) -> Point<'a> {
        let three = FieldElement::new(p1.x.field, 3);
        let slope = (three * (p1.x ^ 2_usize) + self.a) / (p1.y + p1.y);
        let x = (slope ^ 2_usize) - (p1.x + p1.x);
        Point {
            x,
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
                    x: p1.x,
                    y: p1.x.field.zero(),
                }
            }
            return self.point_double(p1);
        }
        let slope = (p2.y - p1.y) / (p2.x - p1.x);
        let x = (slope ^ 2_usize) - p1.x - p2.x;
        Point {
            x,
            y: slope * (p1.x - x),
        }
    }
}
