use crate::perm::Perm;
use num_bigint::BigInt;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Sub, SubAssign};

pub trait ClassLike {
    fn to_repr(&self) -> Perm;
    fn get_class(p: &Perm) -> Self;
    fn get_size(&self) -> BigInt;
}

impl ClassLike for Perm {
    fn to_repr(&self) -> Perm {
        self.clone()
    }
    fn get_class(p: &Perm) -> Self {
        p.clone()
    }
    fn get_size(&self) -> BigInt {
        1.into()
    }
}
pub struct CharEntry {
    pub name: String,
    pub order: BigInt,
    pub table: Vec<BigInt>,
}

impl CharEntry {
    pub fn new<C: ClassLike>(name: String, val: BigInt, g: &[C]) -> CharEntry {
        let mut order = BigInt::from(0);
        let n = g.len();
        for g in g {
            order += g.get_size();
        }
        CharEntry {
            name,
            order,
            table: vec![val; n],
        }
    }
    pub fn from_order(name: String, n: usize, val: BigInt, order: BigInt) -> CharEntry {
        CharEntry {
            name,
            order,
            table: vec![val; n],
        }
    }
    pub fn len(&self) -> usize {
        self.table.len()
    }
    pub fn dim(&self) -> BigInt {
        // Assuming table[0] is the character of e
        self.table[0].clone()
    }
    pub fn inner_prod<C: ClassLike>(&self, other: &Self, g: &[C]) -> BigInt {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.order, other.order);
        let mut sum = BigInt::from(0);
        let order = &self.order;
        for i in 0..self.len() {
            sum += &self.table[i] * &other.table[i] * g[i].get_size();
        }
        assert_eq!(&sum % order, 0.into());
        sum / order
    }
    pub fn tensor_prod(&self, other: &Self) -> Self {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.order, other.order);
        let mut ret = CharEntry::from_order(
            "tensor".to_string(),
            self.len(),
            0.into(),
            self.order.clone(),
        );
        for i in 0..self.len() {
            ret.table[i] = &self.table[i] * &other.table[i];
        }
        ret
    }
    pub fn ext_power_2nd<C: ClassLike + Eq + Hash + Clone>(&self, g: &[C]) -> Self {
        let mut map = HashMap::new();
        let n = g.len();
        for i in 0..n {
            map.insert(g[i].clone(), i);
        }
        let mut ret = CharEntry::from_order(
            "Lambda^2".to_string(),
            self.len(),
            0.into(),
            self.order.clone(),
        );
        for i in 0..self.len() {
            let g2 = g[i].to_repr();
            let g2 = g2.compose(&g2);
            let idx = map[&C::get_class(&g2)];
            ret.table[i] = (&self.table[i] * &self.table[i] - &self.table[idx]) / 2;
        }
        ret
    }
    pub fn ext_power_3rd<C: ClassLike + Eq + Hash + Clone>(&self, g: &[C]) -> Self {
        let mut map = HashMap::new();
        let n = g.len();
        for i in 0..n {
            map.insert(g[i].clone(), i);
        }
        let mut ret = CharEntry::from_order(
            "Lambda^3".to_string(),
            self.len(),
            0.into(),
            self.order.clone(),
        );
        for i in 0..self.len() {
            let g = g[i].to_repr();
            let g2 = g.compose(&g);
            let g3 = g2.compose(&g);
            let idx2 = map[&C::get_class(&g2)];
            let idx3 = map[&C::get_class(&g3)];
            let ch1 = &self.table[i];
            let ch2 = &self.table[idx2];
            let ch3 = &self.table[idx3];
            ret.table[i] = (ch1 * ch1 * ch1 - ch1 * ch2 * 3 + 2 * ch3) / 6;
        }
        ret
    }
    pub fn power(&self, k: BigInt) -> Self {
        let mut ret = CharEntry::from_order(
            "power".to_string(),
            self.len(),
            0.into(),
            self.order.clone(),
        );
        for i in 0..self.len() {
            ret.table[i] = &self.table[i] * &k;
        }
        ret
    }
}

impl Add for &CharEntry {
    type Output = CharEntry;
    fn add(self, other: &CharEntry) -> CharEntry {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.order, other.order);
        let mut ret =
            CharEntry::from_order("+".to_string(), self.len(), 0.into(), self.order.clone());
        for i in 0..self.len() {
            ret.table[i] = &self.table[i] + &other.table[i];
        }
        ret
    }
}

impl AddAssign<&CharEntry> for CharEntry {
    fn add_assign(&mut self, other: &CharEntry) {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.order, other.order);
        for i in 0..self.len() {
            self.table[i] += &other.table[i];
        }
    }
}

impl Sub for &CharEntry {
    type Output = CharEntry;
    fn sub(self, other: &CharEntry) -> CharEntry {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.order, other.order);
        let mut ret =
            CharEntry::from_order("-".to_string(), self.len(), 0.into(), self.order.clone());
        for i in 0..self.len() {
            ret.table[i] = &self.table[i] - &other.table[i];
        }
        ret
    }
}

impl SubAssign<&CharEntry> for CharEntry {
    fn sub_assign(&mut self, other: &CharEntry) {
        assert_eq!(self.len(), other.len());
        assert_eq!(self.order, other.order);
        for i in 0..self.len() {
            self.table[i] -= &other.table[i];
        }
    }
}
