use crate::vc::*;

use groupy::{CurveProjective,CurveAffine};
use paired::{
    bls12_381::*,
    Compress, Engine, PairingCurveAffine,
};
use fff::Field;

use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;

type Idx = usize;
type Val = Xt;

#[derive(Debug)]
pub struct SparseLinFun
{
    pub f: Vec<(Idx, Val)>,
    pub n: usize,
}

impl SparseLinFun {
    pub fn idx(&self, j: usize) -> Idx {
        self.f[j].0
    }

    pub fn val(&self, j: usize) -> Val {
        self.f[j].1
    }

    // (sparse) Inner product
    pub fn ip(&self, v: &Vec<Xt>) -> Xt
    {
        self.f.iter().fold(
            Xt::zero(),
            |sum, &(idx, val)| add_mul!(v[idx], &val, &sum) )
        
    }

    pub fn mul_scal2(&self, g2s: &Vec<G2>) -> G2
    {
        let mut res = G2::zero();
        for (idx, val) in self.f.iter() {
            let mut tmp = g2s[*idx];
            tmp.mul_assign(*val);
            res.add_assign(&tmp);
        }
        res
    }

    pub fn reverse(&self) -> Self
    {
        SparseLinFun {
            f: self.f.iter().map(
                |(idx,val)| (self.n-idx-1, val.clone()))
                .collect(),
            n: self.n
        }
    }

    // multiplies itself to the poly coefficients in v and returns a new vector
    pub fn mul_as_poly(&self, v:&Vec<Xt>) -> Vec<Xt>
    {
        //println!("{} {}", v.len(), self.n);
        let new_poly_size = v.len()+self.n-2+1;
        let mut p:Vec<Xt> = vec![Xt::zero(); new_poly_size];

        for (idx,fval) in self.f.iter() {
            for j in 0..v.len() {
                let mut tmp = v[j];
                tmp.mul_assign(&fval);
                p[idx+j].add_assign(&tmp);
            }
        }
        p
    }


}


// Linear Vector Commitment Trait
pub trait LVC {
    fn commit(pp: &PP, v: &Vec<Xt>) -> Comm;

    fn prove_opn(pp: &PP, v: &Vec<Xt>, f: &SparseLinFun) -> OpnPrf;

    fn vfy_opn(pp: &PP, c: &Comm, f: &SparseLinFun, y: &Xt, prf: &OpnPrf) -> bool;
}

// convert a position to an elementary vector representation
fn vec_fn_from_pos(n: u64, i: u64) -> Vec<Xt>
{
    let mut f: Vec<Xt> = vec![Xt::zero(); n as usize];
    f[i as usize] = Xt::one();
    f
}

// convert a position to an elementary vector representation
fn sparse_lin_fn_from_pos(n: u64, i: u64) -> SparseLinFun
{
    SparseLinFun {
        f: vec![(i as usize, Xt::one())],
        n: n as usize
    }
}

// VC from LVC
impl<T: LVC> VC for T {
    fn commit(pp: &PP, v: &Vec<Xt>) -> Comm
    {
        <T as LVC>::commit(pp, v)
    }

    fn prove_opn(pp: &PP, v: &Vec<Xt>, i: u64) -> OpnPrf
    {
        let f = sparse_lin_fn_from_pos(pp.n, i);
        <T as LVC>::prove_opn(pp, v, &f)
    }

    fn vfy_opn(pp: &PP, c: &Comm, i: u64, x: &Xt, prf: &OpnPrf) -> bool
    {
        let f = sparse_lin_fn_from_pos(pp.n, i);
        <T as LVC>::vfy_opn(pp, c, &f, x, prf)
    }

}
