/* 
Implementation of our submission (Monomial basis) 
*/ 
use crate::vc::*;
use crate::lvc::*;


use groupy::{CurveProjective,CurveAffine};
use paired::{
    bls12_381::*,
    Compress, Engine, PairingCurveAffine,
};

use fff::Field;

use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;

pub struct Mongoose;

// converts a poly in g1
fn poly_in_g1(pp: &PP, v: &Vec<Xt>) -> G1
{
    multi_scalar!(pp.g1sAll, v, 0, v.len())
}


impl LVC for Mongoose {

    fn commit(pp: &PP, v: &Vec<Xt>) -> Comm
    {
        poly_in_g1(pp, v)
    }


    fn prove_opn(pp: &PP, v: &Vec<Xt>, f: &SparseLinFun) -> OpnPrf
    {
        let y = f.ip(v);
        println!("y: {}", y);
        let f_rev = f.reverse();
        let mut mul_poly = f_rev.mul_as_poly(v);

        // this is correct
        println!("f_rev: {:?}, mul_poly: {:?}", f_rev, mul_poly);

        mul_poly[f.n-1].sub_assign(&y);
        println!("mul_poly after subtraction: {:?}", mul_poly);

        // R = first m-1 values of mul_poly
        // H = following m-1 values of mul_poly
        // R_Hat = shift R by one position right
        // just return each of these as multiplication with g1sAll 

        // split
        let H = mul_poly.split_off(f.n);
        let R = mul_poly;
        let mut R_hat = R.clone();
        R_hat.insert(0, Xt::zero()); // multiplies R(X) by X

        println!("H: {:?}, R: {:?}", H, R);


        OpnPrf (
            poly_in_g1(pp, &H),
            poly_in_g1(pp, &R),
            poly_in_g1(pp, &R_hat),
        )

    }

    fn vfy_opn(pp: &PP, c: &Comm, f: &SparseLinFun, y: &Xt, prf: &OpnPrf) -> bool
    {
        let OpnPrf(H, R, R_hat) = prf;

        // degree check
        let mut lhs_deg = Bls12::pairing(R.into_affine(), pp.g2sAll[1].into_affine());
        let rhs_deg =  Bls12::pairing(R_hat.into_affine(), pp.g2);

        let deg_chk_ok =  lhs_deg == rhs_deg;

        println!("deg_chk: {}", deg_chk_ok);

        // for debug
        /*
        let mut tau1_2 = pp.g1sAll[1];
        tau1_2.add_assign(&pp.g1sAll[2]);
        println!("H_eq: {}", tau1_2 == *H);

        let mut tau3_t = Bls12::pairing(pp.g1, pp.g2sAll[3].into_affine());
        let mut tau5_6 = pp.g1sAll[5];
        tau5_6.add_assign(&pp.g1sAll[6]);
        let mut tau5_6_t = Bls12::pairing(tau5_6, pp.g2);

        let mut tau3_5_6 = tau5_6.clone();
        tau3_5_6.add_assign(&pp.g1sAll[3]);
        let mut tau3_5_6_t = Bls12::pairing(tau3_5_6, pp.g2);
*/

        // equation check

        // Cb
        let f_rev = f.reverse();
        let Cb  = f_rev.mul_scal2(&pp.g2sAll);
        //println!("Cb equality {}", Cb == pp.g2sAll[3]); // correct

        /*
        let mut chk_Ca = pp.g1sAll[0];
        chk_Ca.add_assign(&pp.g1sAll[2]);
        chk_Ca.add_assign(&pp.g1sAll[3]);
        println!("Ca equality {}", *c == chk_Ca);
        */

        let CaCb_pair = Bls12::pairing(c.into_affine(), Cb.into_affine());
        //println!("Cabpair equality {}", CaCb_pair == tau3_5_6_t);

        let mut yNMinus1_g1 = pp.g1sAll[f.n-1];
        yNMinus1_g1.mul_assign(*y);
        let mut y_pair = Bls12::pairing(yNMinus1_g1.into_affine(), pp.g2);
        // should be equal to tau^3 in Gt
        //println!("y_pair: {}", y_pair ==  tau3_t);

        let R_pair = Bls12::pairing(R.into_affine(), pp.g2);// should be 0 in Gt
        
        let H_pair = Bls12::pairing(H.into_affine(), pp.g2sAll[f.n]); // should be tau^5 + tau^6 in Gt // and it is
        //println!("H_pair 1: {:?}", H_pair);
        //println!("H_pair 2: {:?}", tau5_6_t);

        //println!("H_pair eq: {}", H_pair ==  tau5_6_t);


        let mut lhs_eq = CaCb_pair; // == tau3_5_6_t
        let y_pair_inv = y_pair.inverse().unwrap();
        lhs_eq.mul_assign(&y_pair_inv); // == tau3_t
        // should be equal to tau^5 + tau^6 in Gt
        //println!("lhs_pair eq: {}", lhs_eq ==  tau5_6_t);

        let mut rhs_eq = R_pair; 
        rhs_eq.mul_assign(&H_pair);

        let eq_chk_ok = lhs_eq == rhs_eq;

        println!("eq_chk: {}", eq_chk_ok);

        deg_chk_ok && eq_chk_ok
    }


}