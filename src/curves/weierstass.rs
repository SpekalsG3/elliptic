use num_bigint::BigUint;
use num_traits::One;
use crate::curves::point::Point;
use crate::field::field_element::FieldElement;
use crate::utils::bit_iter::BitIter;
use crate::utils::s_tonelli::tonelli_shanks;

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

    pub fn evaluate_y(
        &self,
        x: FieldElement<'a>,
    ) -> Option<FieldElement<'a>> {
        tonelli_shanks(
            x.clone() * x.clone() * x.clone()
                + self.a.clone() * x.clone()
                + self.b.clone()
        ).map(|x| x.0)
    }

    pub fn point_double(
        &self,
        p1: Point<'a>,
    ) -> Point<'a> {
        let three = FieldElement::new(p1.x.field, BigUint::from(3_u8));
        let slope = (three.clone() * p1.x.clone() * p1.x.clone() + self.a.clone()) / (p1.y.clone() + p1.y.clone());
        let x = slope.clone() * slope.clone() - p1.x.clone() - p1.x.clone();
        Point {
            x: x.clone(),
            y: (three.clone() * p1.x.clone()) * slope.clone() - (slope.clone() ^ three.value) - p1.y,
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
        let three = FieldElement::new(p1.x.field, BigUint::from(3_u8));
        let slope = (p2.y - p1.y.clone()) / (p2.x.clone() - p1.x.clone());
        let x = slope.clone() * slope.clone() - p1.x.clone() - p2.x.clone();
        Point {
            x: x.clone(),
            y: (p1.x.clone() + p1.x.clone() + p2.x) * slope.clone() - (slope ^ three.value) - p1.y,
        }
    }

    pub fn multi_by_double_and_add(
        &self,
        k: FieldElement<'a>,
        p1: Point<'a>,
    ) -> Point<'a> {
        let mut r0 = p1.clone();
        let mut r1 = p1.clone();

        for digit in k.value.iter_u32_digits() {
            for b in BitIter::from(digit).rev() {
                if b {
                    r0 = self.point_add(r0.clone(), r1.clone());
                }
                r1 = self.point_add(r1.clone(), r1.clone())
            }
        }

        self.point_add(r0, -p1)
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use crate::curves::point::Point;
    use crate::curves::weierstass::WeierstrassCurve;
    use crate::field::field::Field;

    #[test]
    fn test () {
        let field = Field::new(BigUint::from(61_u32));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(9_u8)),
            field.one(),
        );

        let base = {
            let x = field.get(BigUint::from(5_u8));
            Point {
                x: x.clone(),
                y: e.evaluate_y(x).unwrap(),
            }
        };
        assert_eq!(base.y, field.get(BigUint::from(7_u8)));

        let naive_mul = |times: usize| {
            let mut res = base.clone();
            for _ in 1..times {
                res = e.point_add(res, base.clone());
            }
            res
        };

        assert_eq!(e.point_double(base.clone()), Point {
            x: field.get(BigUint::from(26_u8)),
            y: field.get(BigUint::from(50_u8)),
        });
        assert_eq!(naive_mul(3), Point {
            x: field.get(BigUint::from(27_u8)),
            y: field.get(BigUint::from(38_u8)),
        });
        assert_eq!(naive_mul(5), e.multi_by_double_and_add(
            field.get(BigUint::from(5_u8)),
            base.clone(),
        ));
        assert_eq!(e.multi_by_double_and_add(
            field.get(BigUint::from(142_u8)),
            base.clone(),
        ), Point {
            x: field.get(BigUint::from(48_u8)),
            y: field.get(BigUint::from(26_u8)),
        });
    }
}
