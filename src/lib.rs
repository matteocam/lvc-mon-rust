#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(unused_imports)]
#[macro_use] 

mod vc;
mod lvc;

mod mongoose;


macro_rules! gpow_decl {
    ($x:ident, $val:expr) => {
        let mut $x = G1::one();
        $x.mul_assign($val);
    };
}

#[cfg(test)]
mod tests {

    /* test cases:
        - multiscalar DONE
        - committing / proving / verify (correctness) DONE
        - committing / proving / verify (on wrong point) DONE
    */

    use super::*;

    #[test]
    fn setup_vc() {
        use groupy::{CurveProjective,CurveAffine};
        use paired::{
            bls12_381::*,
            Compress, Engine, PairingCurveAffine,
        };
        use fff::Field;

        use vc::*;

        let n = 2;

        let tau_one = Fr::one();
        let mut tau_two = tau_one;
        tau_two.add_assign(&tau_one);

        let pp_one = setup_impl(n, tau_two);

        let tau_np1:Fr = tau_two.pow(&[n+1]); 
        gpow_decl!(g1_np1, tau_np1);

        let exp_val = Bls12::pairing(g1_np1, G2::one());
        assert_eq!(pp_one.gtnp1, exp_val);

    }

    #[test]
    fn grp_ops() {
        use groupy::{CurveProjective,CurveAffine};
        use paired::{
            bls12_381::*,
            Compress, Engine, PairingCurveAffine,
        };
        use fff::Field;

        use vc::*;

        let one_g1 = G1::one();
        let mut two_g1 = one_g1;
        two_g1.add_assign(&one_g1);
        let mut three_g1 = two_g1;
        three_g1.add_assign(&one_g1);

        let mut other_two = three_g1;
        other_two.sub_assign(&one_g1);

        assert_eq!(two_g1, other_two);

        // 1, 2, 1
        let gs = vec![one_g1, two_g1, one_g1];

        // 1, 1, 0
        let ffs = vec![Fr::one(), Fr::one(), Fr::zero()];

        let ms_value:G1 = multi_scalar!(gs, ffs, 0, ffs.len());
        let exp_ms_value = three_g1;
        println!("{}", ms_value);
        assert_eq!(ms_value, exp_ms_value);

        // zero assumption 
        let mut z = G1::zero();
        z.add_assign(&one_g1);
        assert_eq!(z, one_g1);

    }


    #[test]
    fn vc_mongoose() {
        use groupy::{CurveProjective,CurveAffine};
        use paired::{
            bls12_381::*,
            Compress, Engine, PairingCurveAffine,
        };
        use fff::Field;

        use vc::*;
        use mongoose::*;

        let zero = Fr::zero();
        let one = Fr::one();
        // v = [1,0,1,1]
        let v:Vec<Fr> = vec![one, zero, one, one];

        let pp = setup(v.len() as u64);
      
        

        let cm = Mongoose::commit(&pp, &v);

        // test correct value
        {
        let exp_value_idx_0 = one;
        let gd_idx = 0;
        let gd_prf = Mongoose::prove_opn(&pp, &v, gd_idx);
        let gd_result = Mongoose::vfy_opn(&pp, &cm, gd_idx, &exp_value_idx_0, &gd_prf);
        assert!(gd_result); 
        }

        {
        // test wrong value
        let exp_value_idx_1 = zero;
        let bd_idx = 1;
        let bd_prf = Mongoose::prove_opn(&pp, &v, bd_idx);
        let bd_result = Mongoose::vfy_opn(&pp, &cm, bd_idx, &one, &bd_prf);
        assert!(!bd_result); 

        }
    }
}