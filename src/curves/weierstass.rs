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
    order: BigUint,
}

// https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
impl<'a> WeierstrassCurve<'a> {
    pub fn new (
        a: FieldElement<'a>,
        b: FieldElement<'a>,
        order: BigUint,
    ) -> Self {
        assert_eq!(a.field, b.field, "should be in the same field");
        Self {
            a,
            b,
            order,
        }
    }

    pub fn get_base (&self) -> Point {
        let x = self.a.field.get(BigUint::from(5_u8));
        Point {
            x: x.clone(),
            y: self.evaluate_y(x).unwrap(),
            z: self.a.field.one(),
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

    fn point_double(
        &self,
        p1: Point<'a>,
    ) -> Point<'a> {
        if p1.is_infinity() {
            return p1;
        }
        let three = FieldElement::new(p1.x.field, BigUint::from(3_u8));
        let slope = (three.clone() * p1.x.clone() * p1.x.clone() + self.a.clone()) / (p1.y.clone() + p1.y.clone());
        let x = slope.clone() * slope.clone() - p1.x.clone() - p1.x.clone();
        Point {
            z: x.field.one(),
            y: (three.clone() * p1.x.clone()) * slope.clone() - (slope.clone() ^ three.value) - p1.y,
            x,
        }
    }

    pub fn point_add(
        &self,
        p1: Point<'a>,
        p2: Point<'a>,
    ) -> Point<'a> {
        if p1.is_infinity() {
            return p2;
        }
        if p2.is_infinity() {
            return p1;
        }
        if p1.x == p2.x {
            if p1.y == -p2.y {
                return Point::infinity(p1.x.field);
            }
            return self.point_double(p1);
        }
        let three = FieldElement::new(p1.x.field, BigUint::from(3_u8));
        let slope = (p2.y - p1.y.clone()) / (p2.x.clone() - p1.x.clone());
        let x = slope.clone() * slope.clone() - p1.x.clone() - p2.x.clone();
        Point {
            z: x.field.one(),
            y: (p1.x.clone() + p1.x.clone() + p2.x) * slope.clone() - (slope ^ three.value) - p1.y,
            x,
        }
    }

    pub fn project_point_double(
        &self,
        p1: Point<'a>,
    ) -> Point<'a> {
        let two = p1.x.field.get(BigUint::from(2_u8));
        let three = p1.x.field.get(BigUint::from(3_u8));

        let xx = p1.x.clone() * p1.x.clone();
        let zz = p1.z.clone() * p1.z.clone();
        let w = self.a.clone() * zz + three * xx.clone();
        let s = two.clone() * p1.y.clone() * p1.z;
        let ss = s.clone() * s.clone();
        let sss = s.clone() * ss;
        let r = p1.y * s.clone();
        let rr = r.clone() * r.clone();
        let xr = p1.x + r;
        let b = xr.clone() * xr - xx - rr.clone();
        let h = w.clone() * w.clone() - two.clone() * b.clone();

        Point {
            x: h.clone() * s,
            y: w * (b - h) - two * rr,
            z: sss,
        }
    }

    pub fn project_point_add(
        &self,
        p1: Point<'a>,
        p2: Point<'a>,
    ) -> Point<'a> {
        let two = p1.x.field.get(BigUint::from(2_u8));
        let three = p1.x.field.get(BigUint::from(3_u8));

        let u1 = p1.x * p2.z.clone();
        let u2 = p2.x * p1.z.clone();
        let s1 = p1.y * p2.z.clone();
        let s2 = p2.y * p1.z.clone();
        let w = p1.z * p2.z;
        let p = u2.clone() - u1.clone();
        let r = s2.clone() - s1.clone();
        let pp = p.clone() * p.clone();
        let ppp = pp.clone() * p.clone();
        let rr = r.clone() * r.clone();

        Point {
            x: p *(-(u1.clone() + u2.clone()) * pp.clone() + w.clone() * rr.clone()),
            y: (r.clone() * (-two.clone() * w.clone() * rr + three * (u1 + u2) * pp) - ppp.clone() * (s1 + s2)) / two,
            z: w * ppp,
        }
    }

    pub fn double_and_add(
        &self,
        k: BigUint,
        p1: Point<'a>,
    ) -> Point<'a> {
        let infinity = Point::infinity(p1.x.field);
        let mut r0 = infinity.clone();
        let mut r1 = p1;

        for digit in k.iter_u32_digits() {
            for b in BitIter::from(digit).rev() {
                if b {
                    r0 = self.point_add(r0.clone(), r1.clone());
                } else {
                    r0 = self.point_add(r0, infinity.clone()); // optional, against cache-timing attack
                }
                r1 = self.point_add(r1.clone(), r1.clone());
            }
        }

        r0
    }

    pub fn montgomery_ladder(
        &self,
        k: BigUint,
        p1: Point<'a>,
    ) -> Point<'a> {
        // find `n` where `curve order < 2 ^ n`
        let size = u32::BITS as usize;
        let mut bits = 0;
        {
            let mut order_iter = self.order.iter_u32_digits().rev();
            if let Some(digit) = order_iter.next() {
                let bit_iter = BitIter::from(digit).skip_while(|b| b == &false);
                bits += bit_iter.count();
            }
            bits += order_iter.count() * size;
        }

        // convert k to keep `nth` bit always at one
        // this lets us avoid setting r0 to point of infinity and handling it during operations
        // by implicitly performing one iteration
        let k = {
            let n_2 = BigUint::from(1_u8) << bits;
            let d = if k > n_2 {
                (k - n_2.clone()) % self.order.clone()
            } else {
                self.order.clone() - (n_2.clone() - k) % self.order.clone()
            };

            n_2 + d
        };

        let mut r0 = p1.clone();
        let mut r1 = self.project_point_double(p1.clone());

        for digit in k.iter_u32_digits().rev() {
            let count = bits % size;
            bits -= count;

            for b in BitIter::from(digit).skip(size - count) {
                if b {
                    r0 = self.project_point_add(r0.clone(), r1.clone());
                    r1 = self.project_point_double(r1);
                } else {
                    r1 = self.project_point_add(r0.clone(), r1.clone());
                    r0 = self.project_point_double(r0);
                }
            }
        }

        r0
    }

    pub fn primitive_root (&self) -> Point<'_> {
        let y = -self.a.field.get(BigUint::from(2_u8));

        // let x = (-self.get(BigUint::from(2_u8))) ^ self.order.clone();
        let x = tonelli_shanks(-self.a.field.get(BigUint::from(3_u8)));
        if let Some(x) = x {
            // y = ABS(x*y - y,p);
            let y = x.0 * y.clone() - y;

            // return newfpoint(0,y);
            return Point {
                x: self.a.field.zero(),
                y,
                z: self.a.field.one(),
            };
        }

        let yy = y.clone() * y.clone();
        let x = yy - y.clone();
        let x = x + self.a.field.one();

        // x = modsquareroot(ABS(y*y - y + 1,p),p);
        let x = tonelli_shanks(x).expect("should return proper sqrt");

        // return newfpoint(x,ABS(-y,p));
        Point { x: x.0, y: -y, z: self.a.field.one() }
    }

    // order of a subgroup generated by provided point
    pub fn find_order(
        &self,
        mut p: Point<'a>,
    ) -> BigUint {
        if p.is_infinity() {
            return BigUint::one();
        }

        let mut n = BigUint::from(2_u8);
        let mut t = self.project_point_double(p.clone());
        while !t.is_infinity() {
            t = self.project_point_add(t, p.clone());
            n += BigUint::one();
        }

        n
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::One;
    use crate::curves::point::Point;
    use crate::curves::weierstass::WeierstrassCurve;
    use crate::field::field::Field;

    #[test]
    fn get_order() {
        let field = Field::new(BigUint::from(631_u16));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(9_u8)),
            field.one(),
            BigUint::from(73_u8),
            field.get(BigUint::from(30_u8)),
            field.get(BigUint::from(34_u8)),
            BigUint::from(5_u8),
        );


        {
            let p = Point {
                x: field.get(BigUint::from(36_u8)),
                y: field.get(BigUint::from(60_u8)),
                z: field.one(),
            };
            assert_eq!(e.find_order(p), (5_u8).into());
        }
        {
            let q = Point {
                x: field.get(BigUint::from(121_u8)),
                y: field.get(BigUint::from(387_u16)),
                z: field.one(),
            };
            assert_eq!(e.find_order(q), (5_u8).into());
        }
        {
            let s = Point {
                x: field.get(BigUint::from(0_u8)),
                y: field.get(BigUint::from(36_u16)),
                z: field.one(),
            };
            assert_eq!(e.find_order(s), (130_u8).into());
        }
    }

    #[test]
    fn primitive_nth_root () {
        let p: u8 = 61;
        let curve_order: usize = 73;

        let field = Field::new(BigUint::from(p));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(9_u8)),
            field.one(),
            BigUint::from(curve_order),
        );

        let z = e.primitive_root();

        let first = e.montgomery_ladder(BigUint::one(), z.clone());
        let mut n = 2;
        while n < curve_order+10 {
            let r = e.montgomery_ladder(BigUint::from(n), z.clone());
            if r == first {
                break
            }
            n += 1;
        }
        assert_eq!(n - 1, curve_order);
    }

    #[test]
    fn test () {
        let field = Field::new(BigUint::from(61_u8));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(9_u8)),
            field.one(),
            BigUint::from(73_u8),
        );

        let base = e.get_base();
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
            z: field.one(),
        });
        assert_eq!(naive_mul(3), Point {
            x: field.get(BigUint::from(27_u8)),
            y: field.get(BigUint::from(38_u8)),
            z: field.one(),
        });
        assert_eq!(naive_mul(5), e.double_and_add(
            BigUint::from(5_u8),
            base.clone(),
        ));
        assert_eq!(e.double_and_add(
            BigUint::from(142_u8),
            base.clone(),
        ), Point {
            x: field.get(BigUint::from(48_u8)),
            y: field.get(BigUint::from(26_u8)),
            z: field.one(),
        });
        {
            let p = e.project_point_double(base.clone()); // 2
            assert_eq!(p.x.clone() / p.z.clone(), field.get(BigUint::from(26_u8)));
            assert_eq!(p.y.clone() / p.z.clone(), field.get(BigUint::from(50_u8)));

            let p = e.project_point_add(p, base.clone()); // 3
            assert_eq!(p.x.clone() / p.z.clone(), field.get(BigUint::from(27_u8)));
            assert_eq!(p.y.clone() / p.z.clone(), field.get(BigUint::from(38_u8)));

            let p = e.project_point_add( // 5
                e.project_point_double(p), // 6
                -base.clone(), // -1
            );
            assert_eq!(p.x.clone() / p.z.clone(), field.get(BigUint::from(30_u8)));
            assert_eq!(p.y.clone() / p.z.clone(), field.get(BigUint::from(59_u8)));
        }

        {
            let p = e.montgomery_ladder(BigUint::from(5_u8), base.clone());
            assert_eq!(p.x / p.z.clone(), field.get(BigUint::from(30_u8)));
            assert_eq!(p.y / p.z.clone(), field.get(BigUint::from(59_u8)));

            let p = e.montgomery_ladder(BigUint::from(142_u8), base.clone());
            assert_eq!(p.x / p.z.clone(), field.get(BigUint::from(48_u8)));
            assert_eq!(p.y / p.z.clone(), field.get(BigUint::from(26_u8)));
        }
    }
}
