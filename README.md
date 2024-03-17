# Project

This project uses same `FieldElement` and `Field` defined in `zk-stark-tutor`,
except one change - in elliptic curve arithmetics,
both identity element and zero element is so-called point of infinity,
which is why `Neg` implementation is removed to avoid confusion.

# Sources

- https://curves.xargs.org/
- https://www.hyperelliptic.org/EFD/
- https://martin.kleppmann.com/papers/curve25519.pdf
- https://eprint.iacr.org/2020/437.pdf
- https://electricdusk.com/files/curve13318_indocrypt19-20191218.pdf

- https://static1.squarespace.com/static/5fdbb09f31d71c1227082339/t/5ff394720493bd28278889c6/1609798774687/PairingsForBeginners.pdf
