use num_bigint::BigUint;
use num_traits::One;
use rand::{RngCore, thread_rng};
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

    pub fn get_base (&self) -> Point {
        let x = self.a.field.get(BigUint::from(5_u8));
        Point {
            x: x.clone(),
            y: self.evaluate_y(x).unwrap().0,
            z: self.a.field.one(),
        }
    }

    pub fn evaluate_y(
        &self,
        x: FieldElement<'a>,
    ) -> Option<(FieldElement<'a>, FieldElement<'a>)> {
        tonelli_shanks(
            x.clone() * x.clone() * x.clone()
                + self.a.clone() * x.clone()
                + self.b.clone()
        )
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
        let three = BigUint::from(3_u8);
        let slope = (p2.y - p1.y.clone()) / (p2.x.clone() - p1.x.clone());
        let x = slope.clone() * slope.clone() - p1.x.clone() - p2.x.clone();
        Point {
            z: x.field.one(),
            y: (p1.x.clone() + p1.x.clone() + p2.x) * slope.clone() - (slope ^ three) - p1.y,
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
        let order = self.find_order(p1.clone());

        let size = u32::BITS as usize;
        let mut bits = 0;
        {
            let mut order_iter = order.iter_u32_digits().rev();
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
                (k - n_2.clone()) % order.clone()
            } else {
                order.clone() - (n_2.clone() - k) % order.clone()
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

    // order of a subgroup generated by provided point
    pub fn find_order(
        &self,
        p: Point<'a>,
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

    // evaluate addition of
    // chord line going through `P` and `Q` (or tangent line if `P==Q`)
    // and vertical line going through `P+Q` and `-(P+Q)`
    // in the point `S`
    pub fn eval_chord_tangent(
        &self,
        p: Point<'a>,
        q: Point<'a>,
        s: Point<'a>,
    ) -> Option<FieldElement<'a>> {
        if p.is_infinity() || q.is_infinity() {
            return None;
        }

        let result;
        if p == -q.clone() { // when slope is infinity
            result = s.x - p.x;
        } else {
            let slope = if p.x == q.x {
                let three = self.a.field.get(BigUint::from(3_u8));
                (three * p.x.clone() * p.x.clone() + self.a.clone()) / (p.y.clone() + p.y.clone())
            } else {
                (q.y - p.y.clone()) / (q.x.clone() - p.x.clone())
            };
            let num = s.y - p.y - slope.clone() * (s.x.clone() - p.x.clone());
            let din = s.x + p.x + q.x - slope.clone() * slope;
            result = num / din;
        }

        Some(result)
    }

    // evaluation of Miller function `f_{m,P}` in point `Q` (rational function satisfying `div(f_{m,P}) = m[P] − [mP] − (m − 1)[O]`)
    fn miller(
        &self,
        p: Point<'a>,
        q: Point<'a>,
        m: BigUint,
    ) -> Option<FieldElement<'a>> {
        let size = u32::BITS as usize;
        let mut bits = 0;
        {
            let mut order_iter = m.iter_u32_digits().rev();
            if let Some(digit) = order_iter.next() {
                let bit_iter = BitIter::from(digit).skip_while(|b| b == &false);
                bits += bit_iter.count();
            }
            bits += order_iter.count() * size;
        }
        let array = m
            .iter_u32_digits()
            .flat_map(|digit| {
                let count = bits % size;
                bits -= count;
                BitIter::from(digit).skip(size - count)
            });

        let mut t = p.clone();
        let mut f = self.a.field.one();

        // l_{[m]*T,T} / v_{[m+1]*T} where
        // `l` is the line through `T` and `[m]T` (or tangent if `T=[m]T`),
        // and `v` is the vertical through `[m+1]T` and `-[m+1]T`
        // has the divisor `f_{m+1,T} − f_{m,T}`, thus we use that to obtain `f_{m+1,T} = f_{m,T} * g_{m,T}`
        // and similarly `l_{[m]T,[m]T} / v_{[2m]T}` has divisor for `f_{2m,T} - f_{m,T}^2`
        // so for any `m` we can quickly obtain `f_{m+1,T}` and `f_{2m,T}`
        // this is where the double-and-add-like Miller algorithm comes from

        for byte in array.skip(1) {
            let temp = self.eval_chord_tangent(t.clone(), t.clone(), q.clone())?; // eval `l_{[m]*T,T} / v_{[m+1]*T}` in `Q`
            f = f.clone() * f; // `f_{m,P} ^ 2`
            f = f * temp;
            t = self.point_double(t); // t = 2t // todo: project_point_double

            if byte {
                // g_{t,P}(Q) = l is the line through t and a, l' through t+a and -(t+P)
                let temp = self.eval_chord_tangent(t.clone(), p.clone(), q.clone())?; // eval `l_{[m]*T,T} / v_{[m+1]*T}` in `Q`
                f = f * temp; // f = f * g_{t,P}(Q)
                t = self.point_add(t, p.clone()); // t = t + P // todo: project_point_add
            }
        }

        Some(f)
    }

    pub fn random_point(
        &self,
    ) -> Option<Point<'a>> {
        let mut thread_rng = thread_rng();
        let mut bytes = vec![0; self.a.field.order.bits() as usize];

        for _ in 0..100 {
            thread_rng.fill_bytes(&mut bytes);
            let x = self.a.field.sample(&bytes);

            if let Some(y) = self.evaluate_y(x.clone()) {
                let s = Point {
                    x,
                    y: y.0,
                    z: self.a.field.one(),
                };

                return Some(s);
            }
        }
        None
    }

    fn phi(
        &self,
        a: Point<'a>,
    ) -> Point<'a> {
        if a.is_infinity() {
            return a;
        }

        // FPOINT * temp = primitroot(p);
        let temp = self.a.field.generator().expect("should have generator");

        // fmulti(a->x,temp,p,result->x);
        let x = a.x * temp;
        // assign(result->y, a->y);
        let y = a.y;

        Point {
            x,
            y,
            z: self.a.field.one(),
        }
    }

    pub fn weilpairing(
        &self,
        m: BigUint,
        p: Point<'a>,
        q: Point<'a>,
    ) -> Option<FieldElement<'a>> {
        // { // we calculate m using the same function, is there a point to do this check?
        //     // Let `P,Q in E[m]` be points of order `m` in `E`
        //     let m_p = self.find_order(p.clone());
        //     let m_q = self.find_order(q.clone());
        //     if m_p != m_q || m_p != m {
        //         return None;
        //     }
        // }

        let mut s = self.random_point().expect("failed to gen point");
        {
            let tries = 100;
            let mut i = 0;
            loop {
                // from definition `S` is any point on the curve not in `{O, P, −Q, P − Q}`
                // but it shows there's also some another property to points `P`,`Q` and `S`
                // because pairing always works if order of S != m
                // but only sometimes works when ord(S) == m (for example, S=(339;499;1))
                // and other times results in division by zero
                // specifically, in evaluation of `l_{[m]P,[m]P} / v_{[2m]P}` in point `Q+S`
                // idk why
                if !s.is_infinity()
                    && s != p
                    && s != -q.clone()
                    && s != self.point_add(p.clone(), -q.clone())
                    && self.find_order(s.clone()) != m
                {
                    break
                }
                s = self.random_point().expect("failed to gen point");
                i += 1;
            }
            if i == tries {
                panic!("failed to generate proper s for provided M");
                // return None;
            }
        }

        let mp_q_s  = self.miller(p.clone(), self.point_add(q.clone(),  s.clone()), m.clone())?;
        let mp_s    = self.miller(p.clone(),                            s.clone() , m.clone())?;
        let mq_p_ns = self.miller(q.clone(), self.point_add(p.clone(), -s.clone()), m.clone())?;
        let mq_ns   = self.miller(q.clone(),                           -s.clone() , m.clone())?;

        // it is `e_m(P,S)'` as inverse of weil pairing `e_m(P,S)` and, in fact, they are equal,
        // but it has all the same properties so application-wise we can use any of them
        Some((mp_q_s / mp_s) / (mq_p_ns / mq_ns))
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use crate::curves::point::Point;
    use crate::curves::weierstass::WeierstrassCurve;
    use crate::field::field::Field;

    #[test]
    fn eval_divisor() {
        let field = Field::new(BigUint::from(61_u8));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(9_u8)),
            field.one(),
        );

        {
            let a = e.get_base();
            let b = e.double_and_add(BigUint::from(25_u32), a.clone());
            let c = e.double_and_add(BigUint::from(53_u32), a.clone());
            assert_eq!(e.eval_chord_tangent(a, b, c), Some(field.get(BigUint::from(12_u8))));
        }
        {
            let a = e.get_base();
            let b = e.double_and_add(BigUint::from(678234_u32), a.clone());
            let c = e.double_and_add(BigUint::from(346857_u32), a.clone());
            assert_eq!(e.eval_chord_tangent(a, b, c), Some(field.get(BigUint::from(20_u8))));
        }
        {
            let a = e.get_base();
            let b = e.double_and_add(BigUint::from(111_u32), a.clone());
            let c = e.double_and_add(BigUint::from(999_u32), a.clone());
            assert_eq!(e.eval_chord_tangent(a, b, c), Some(field.get(BigUint::from(26_u8))));
        }
    }

    #[test]
    fn miller() {
        let field = Field::new(BigUint::from(631_u16));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(30_u8)),
            field.get(BigUint::from(34_u8)),
        );

        let p = Point {
            x: field.get(BigUint::from(36_u8)),
            y: field.get(BigUint::from(60_u8)),
            z: field.one(),
        };
        let q = Point {
            x: field.get(BigUint::from(121_u8)),
            y: field.get(BigUint::from(387_u16)),
            z: field.one(),
        };
        let s = Point {
            x: field.get(BigUint::from(0_u8)),
            y: field.get(BigUint::from(36_u16)),
            z: field.one(),
        };
        let m = BigUint::from(5_u8);
        assert_eq!(
            e.miller(p.clone(), e.point_add(q.clone(), s.clone()), m.clone()),
            Some(field.get(BigUint::from(103_u8))),
        );
        assert_eq!(
            e.miller(p.clone(), s.clone(), m.clone()),
            Some(field.get(BigUint::from(219_u8))),
        );
        assert_eq!(
            e.miller(q.clone(), e.point_add(p, -s.clone()), m.clone()),
            Some(field.get(BigUint::from(284_u16))),
        );
        assert_eq!(
            e.miller(q.clone(), -s.clone(), m.clone()),
            Some(field.get(BigUint::from(204_u8))),
        );
    }

    #[test]
    fn get_order() {
        let field = Field::new(BigUint::from(631_u16));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(30_u8)),
            field.get(BigUint::from(34_u8)),
        );

        {
            let p = Point {
                x: field.get(BigUint::from(36_u8)),
                y: field.get(BigUint::from(60_u8)),
                z: field.one(),
            };
            let m = e.find_order(p.clone());
            assert_eq!(m, (5_u8).into());
            assert_eq!(e.double_and_add(m, p), Point::infinity(&field));
        }
        {
            let q = Point {
                x: field.get(BigUint::from(121_u8)),
                y: field.get(BigUint::from(387_u16)),
                z: field.one(),
            };
            let m = e.find_order(q.clone());
            assert_eq!(m, (5_u8).into());
            assert_eq!(e.double_and_add(m, q), Point::infinity(&field));
        }
        {
            let s = Point {
                x: field.get(BigUint::from(0_u8)),
                y: field.get(BigUint::from(36_u16)),
                z: field.one(),
            };
            let m = e.find_order(s.clone());
            assert_eq!(m, (130_u8).into());
            assert_eq!(e.double_and_add(m, s), Point::infinity(&field));
        }
    }

    #[test]
    pub fn weilpairing() {
        let field = Field::new(BigUint::from(631_u16));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(30_u8)),
            field.get(BigUint::from(34_u8)),
        );

        let p = Point {
            x: field.get(BigUint::from(36_u16)),
            y: field.get(BigUint::from(60_u16)),
            z: field.one(),
        };
        let q = Point {
            x: field.get(BigUint::from(121_u16)),
            y: field.get(BigUint::from(387_u16)),
            z: field.one(),
        };
        let m = e.find_order(p.clone());
        // any `m` will be co-prime to `char(F)` which is prime, no need to check that
        // todo
        //  check `Q` is not a multiple of point `P` to ensure that functions
        //  used during evaluation of miller function are not vanished
        //  otherwise this will produce degenerative pairing (resulting in `Some(1)`)
        assert_eq!(m, e.find_order(q.clone()));

        let pairing = e.weilpairing(m.clone(), p.clone(), q.clone());

        // properties
        { // non-degenerate
            // assert_eq!(pairing, Some(field.get(BigUint::from(242_u8))));
            assert!(pairing.is_some());
            assert_ne!(pairing, Some(field.one()));
        }
        let pairing = pairing.unwrap();
        { // mth root of unity
            assert_eq!(pairing.clone() ^ m.clone(), field.one());
        }
        { // bilinear
            let a = BigUint::from(7_u8);
            assert_eq!(
                e.weilpairing(m.clone(), e.double_and_add(a.clone(), p.clone()), q.clone()),
                Some(pairing.clone() ^ a.clone()),
            );
            assert_eq!(
                e.weilpairing(m.clone(), p.clone(), e.double_and_add(a.clone(), q.clone())),
                Some(pairing.clone() ^ a.clone()),
            );
        }
        {// alternating
            let a = BigUint::from(7_u8);
            assert_eq!(
                e.weilpairing(m.clone(), p.clone(), p.clone()),
                Some(field.one()),
            );
            let p2 = e.double_and_add(a.clone(), p.clone());
            assert_eq!(
                e.weilpairing(m.clone(), p2.clone(), p2.clone()),
                Some(field.one()),
            );
        }
    }

    #[test]
    fn adds_n_mults() {
        let field = Field::new(BigUint::from(61_u8));
        let e = WeierstrassCurve::new(
            field.get(BigUint::from(9_u8)),
            field.one(),
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
