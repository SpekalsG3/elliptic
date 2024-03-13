use crate::field::field_element::FieldElement;

/// Tonelli-Shanks algorithm
/// Find quadratic residue `n` for given `x`, such that `x^2 = n mod p`, where p is prime
/// According to Euler's criterion, in such field root exists iff `n^{(p-1)/2} = 1 mod p`
pub fn tonelli_shanks(x: FieldElement<'_>) -> Option<(FieldElement<'_>, FieldElement<'_>)> {
    let one = x.field.one();
    let one_inv = -one; // p - 1

    let p = x.field.order;
    let mut q = p - 1;
    let mut ss: u128 = 0;
    let mut z = FieldElement::new(x.field, 2);

    if x ^ ((p - 1) / 2) != one {
        return None;
    }

    while (q & 1) == 0 {
        ss += 1;
        q >>= 1;
    }

    if ss == 1 {
        let r1 = x ^ ((p + 1) / 4); // +1 ???
        return Some((r1, -r1));
    }

    while z ^ ((p - 1) / 2) != one_inv {
        z = z + one;
    }

    let mut c = z ^ q;
    let mut r = x ^ ((q + 1) / 2);
    let mut t = x ^ q;
    let mut m = ss;

    loop {
        if t == one {
            return Some((r, -r));
        }
        let mut i = 0;
        let mut zz = t;
        while zz != one && i < (m - 1) {
            zz = zz * zz;
            i += 1;
        }
        let mut b = c;
        let mut e = m - i - 1;
        while e > 0 {
            b = b * b;
            e -= 1;
        }
        c = b * b;
        t = t * c;
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
            let field = Field::new(p);
            assert_eq!(
                tonelli_shanks(FieldElement::new(&field, n)),
                r.map(|(r1, r2)| (FieldElement::new(&field, r1), FieldElement::new(&field, r2)))
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
