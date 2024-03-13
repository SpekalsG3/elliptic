# Project

This project uses same `FieldElement` and `Field` defined in `zk-stark-tutor`,
except one change - in elliptic curve arithmetics,
both identity element and zero element is so-called point of infinity,
which is why `Neg` implementation is removed to avoid confusion.
