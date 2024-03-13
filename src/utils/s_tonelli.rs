use num_bigint::BigUint;
use num_traits::Zero;
use crate::field::field_element::FieldElement;

/// Tonelli-Shanks algorithm
/// Find quadratic residue `n` for given `x`, such that `x^2 = n mod p`, where p is prime
/// According to Euler's criterion, in such field root exists iff `n^{(p-1)/2} = 1 mod p`
pub fn tonelli_shanks(x: FieldElement<'_>) -> Option<(FieldElement<'_>, FieldElement<'_>)> {
    let big_one = BigUint::from(1_u8);
    let big_two = BigUint::from(2_u8);

    let one = x.field.one();
    let one_inv = -one.clone(); // p - 1

    let p = x.field.order.clone();
    let mut q = p.clone() - big_one.clone();
    let mut ss: u128 = 0;
    let mut z = FieldElement::new(x.field, big_two.clone());

    if x.clone() ^ ((p.clone() - big_one.clone()) / big_two.clone()) != one {
        return None;
    }

    while (q.clone() & big_one.clone()).is_zero() {
        ss += 1;
        q >>= 1;
    }

    if ss == 1 {
        let r1 = x.clone() ^ ((p.clone() + big_one.clone()) / BigUint::from(4_u8)); // +1 ???
        return Some((r1.clone(), -r1));
    }

    while z.clone() ^ ((p.clone() - big_one.clone()) / big_two.clone()) != one_inv {
        z = z.clone() + one.clone();
    }

    let mut c = z ^ q.clone();
    let mut r = x.clone() ^ ((q.clone() + big_one) / big_two);
    let mut t = x ^ q;
    let mut m = ss;

    loop {
        if t == one {
            return Some((r.clone(), -r));
        }
        let mut i = 0;
        let mut zz = t.clone();
        while zz != one && i < (m - 1) {
            zz = zz.clone() * zz;
            i += 1;
        }
        let mut b = c.clone();
        let mut e = m - i - 1;
        while e > 0 {
            b = b.clone() * b;
            e -= 1;
        }
        c = b.clone() * b.clone();
        t = t * c.clone();
        r = r * b;
        m = i;
    }
}

#[cfg(test)]
mod tests {
    use crate::field::field::Field;
    use crate::field::field_element::FieldElement;
    use crate::utils::s_tonelli::tonelli_shanks;

    #[test]
    fn test () {
        fn run (p: u128, n: u128, r: Option<(u128, u128)>) {
            let field = Field::new(p.into());
            assert_eq!(
                tonelli_shanks(FieldElement::new(&field, n.into())),
                r.map(|(r1, r2)| (FieldElement::new(&field, r1.into()), FieldElement::new(&field, r2.into())))
            )
        }

        run(13, 10, Some((7, 6)));
        run(101, 56, Some((37, 64)));
        run(10009, 1030, Some((1632, 8377)));
        run(100049, 44402, Some((30468, 69581)));
        run(1000000009, 665820697, Some((378633312, 621366697)));
        run(1000000000039, 881398088036, Some((791399408049, 208600591990)));
    }
}
