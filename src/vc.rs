use groupy::{CurveProjective,CurveAffine};
use paired::{
    bls12_381::*,
    Compress, Engine, PairingCurveAffine,
};
use fff::Field;

use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;

// pub type for values we are committing to
pub type Xt = Fr;
pub type Gt = <Bls12 as Engine>::Fqk;

// Simplifying assumption (valid for Pointproofs and our conpub struction):
//  an opening proof is either a group element or a tuple of group elements
pub struct OpnPrf(pub G1, pub G1, pub G1);

impl OpnPrf {
    pub fn just_one(g: G1) -> OpnPrf
    {
        OpnPrf(g, G1::zero(), G1::zero())
    }
}

pub struct PP {
    pub g1taus: Vec<G1>,
    pub g1sAll: Vec<G1>,
    pub lambdas: Vec<G1>,
    pub g2taus: Vec<G2>,
    pub g2sAll: Vec<G2>,
    pub g1: G1,
    pub g2: G2,
    pub gtnp1: Gt,
    pub n: u64,
}

// auxiliary trait to easily access zeros of group without passing type to macros
pub trait GrpZeroAux
{
    fn grp_zero(&self) -> Self;
}

impl GrpZeroAux for G1
{
    fn grp_zero(&self) -> Self
    {
        G1::zero()
    }
}

impl GrpZeroAux for G2
{
    fn grp_zero(&self) -> Self
    {
        G2::zero()
    }
}

/* 
Why the macros below?
They are a way to compensate for fff's interface. 
All arithmetic operations seem to have side-effects there. 
These macros provide functional versions of them
*/

// returns (sum+a*b)
macro_rules! add_mul {
    ($a:expr, $b:expr , $sum:expr) => {{
        let mut tmp = $a;
        tmp.mul_assign($b);
        tmp.add_assign($sum);
        tmp
    }}
}

// performs sum_i gs[i]*xs[i]
macro_rules! multi_scalar {
    ($gs:expr, $xs:ident, $start:expr, $end: expr) => {{
        let mut res = $gs[0].grp_zero();
        for i in $start..$end {
            let mut tmp = $gs[i];
            tmp.mul_assign($xs[i]);
            res.add_assign(&tmp);
        }
        res
    }};
}

// performs sum_i gs[i+offset]*xs[i]
macro_rules! multi_scalar_offset {
    ($gs:expr, $xs:ident, $start:expr, $end: expr, $offset: expr) => {{
        let mut res = $gs[0].grp_zero();
        for i in $start..$end {
            let mut tmp = $gs[i+$offset];
            tmp.mul_assign($xs[i]);
            res.add_assign(&tmp);
        }
        res
    }};
}


pub fn to_lagrange_basis_grp(gs: Vec<G1>) -> Vec<G1>
{
    unimplemented!();
}

// setup implementation
pub fn setup_impl(n: u64, tau: Fr) -> PP
{
    
    let g1 = G1::one(); 
    let g2 = G2::one();

    let mut g1pow = g1;
    let mut g2pow = g2;

    let mut g1s = Vec::<G1>::new();
    let mut g2s = Vec::<G2>::new();

    // from zero
    let mut g1sAll = Vec::<G1>::new();
    let mut g2sAll = Vec::<G2>::new();
    g1sAll.push(g1);
    g2sAll.push(g2);



    // XXX: There may be a problem with where you are starting from index-wise in the g1s (but the problem might just be in Pointproofs)
    // anyways for now: g1^tau == g1s[0] =/= g1     

    // dummy init value
    let mut gtnp1 = Bls12::pairing(g1, g2); 
    
    // We make implementation easier by setting N+1th power to zero

    for i in 1..2*n {
        g1pow.mul_assign(tau);
        g2pow.mul_assign(tau);

        g1sAll.push(g1pow);
        g2sAll.push(g2pow);
        g2s.push(g2pow);

        if i == n+1 {
            g1s.push(G1::zero());
            //let tau_np1 = tau.pow(&[n+1]); 
            // gtnp1 := [tau^(n+1)]
            gtnp1 = Bls12::pairing(g1pow, g2);
        } else {
            g1s.push(g1pow);
        }
    }

    //g1sAll.extend(g1s.clone());
    //g2sAll.extend(g2s.clone());

    // TODO: add rest
    PP { g1sAll: g1sAll,
         g2sAll: g2sAll,
         g1taus: g1s.clone(),
         g2taus: g2s,
         g1: g1, 
         g2: g2,
         gtnp1: gtnp1,
         lambdas: g1s, //lambdas: to_lagrange_basis_grp(g1s),
         n : n }
}

// Common setup
pub fn setup(n: u64) -> PP
{
    // FIXME: Change this to a different RNG
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    // sample tau
    let tau = Fr::random(&mut rng);

    setup_impl(n, tau)
}

pub type Comm = G1;


pub trait VC {
    fn commit(pp: &PP, v: &Vec<Xt>) -> Comm;

    fn prove_opn(pp: &PP, v: &Vec<Xt>, i: u64) -> OpnPrf;

    fn vfy_opn(pp: &PP, c: &Comm, i: u64, x: &Xt, prf: &OpnPrf) -> bool;
}

