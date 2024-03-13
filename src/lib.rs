
pub mod field;
pub mod utils;
pub mod curves;


#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use std::str::FromStr;

    #[test]
    fn bigint () {
        let x = BigUint::from_str("57896044618658097711785492504343953926634992332820282019728792003956564819949").unwrap();
        println!("{:?}", x.iter_u32_digits().collect::<Vec<_>>());
    }
}
