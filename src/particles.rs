
use std::fmt;
use rand::{distributions::{Distribution, Standard},Rng};

#[derive(Clone)]
pub struct Link{
    target_id:u32,
    coupling:f64,
}

impl Link{
    pub fn new(target_id:u32, coupling:f64)->Self{
        return Link {
             target_id,
             coupling};
    }

    pub fn new_from(target_id:u32, other_link:&Link)->Self{
        return Link{
            target_id: target_id,
            ..*other_link
        }
    }
    pub fn target_id(&self)->u32{return self.target_id}
    pub fn coupling(&self)->f64{return self.coupling}
}

impl PartialEq for Link{
    fn eq(&self, other: &Self) -> bool {
        return self.coupling == other.coupling;
    }
}

#[derive(Clone)]
pub enum Spin{
    Up,
    Down,
}

impl Spin{
    pub fn flip(&self)->Spin{
        match self{
            Spin::Up => return Spin::Down,
            Spin::Down => return Spin::Up,
        }
    }

    pub fn value(&self)->i64{
        match self{
            Spin::Up => return 1,
            Spin::Down => return -1,
        }
    }
}

impl fmt::Display for Spin {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,"{}",self.value().to_string())
    }
}

impl Distribution<Spin> for Standard{
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Spin {
        match rng.gen_bool(0.5){
            true => return Spin::Up,
            false => return Spin::Down,
        }
    }
}

impl From<&Spin> for f64{
    fn from(spin: &Spin) -> Self {
        match spin{
            Spin::Up => return 1.0,
            Spin::Down => return -1.0,
        }
    } 
}




#[derive(Clone)]
pub struct Particle{
    id: u32,
    spin:Spin,
    links: Vec<Link>,
}

impl Particle{
    pub fn new(id:u32, spin:Spin)->Self{
        Particle{
            id,
            spin,
            links:Vec::new(),
        } 
    }

    pub fn remove_neighbour(&mut self, neighbour_id: u32)->&mut Self{
        //Retains all elements in self.neighbours not equal to neighbour_idx
        //  -> all equal get removed
        self.links.retain(|link| link.target_id != neighbour_id);
        return self
    }
    pub fn link_neighbour(&mut self, link: Link)->&mut Self{
        //link neighbour by pushing its idx in grid to self.neighbours vec
        if self.id == link.target_id{
            panic!("You tried to link a particle to itself");
        }
        // if link already exists in particle
        if self.find_link(link.target_id).is_some(){
            return self
        }
        self.links.push(link);
        return self
    }
    pub fn link_neighbours(&mut self, link_vec: &Vec<Link>){
        for link in link_vec.iter(){
            self.link_neighbour(link.clone());
        }
    }

    pub fn id(&self)->u32{return self.id}
    pub fn neighbours(&self)->&Vec<Link>{return &self.links}
    pub fn find_link(&self, link_idx:u32)->Option<Link>{
        // tries to find a existing Link. Returns Some(copy) of link if found.
        match self.links.iter().find(|&link|link.target_id == link_idx){
            None => return None,
            Some(link) => return Some(link.clone()),
        }
        
    }
    pub fn spin(&self)->&Spin{return &self.spin}
    pub fn flip_spin(&mut self){self.spin = self.spin.flip();}
}

impl fmt::Display for Particle{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //Neighbours: 1,2,3,4
        let mut output = String::new();
        for neighbour in self.links.iter(){
            output.push_str(&format!("
                (idx: {}, coupling: {})",
                    &neighbour.target_id.to_string(),
                    &neighbour.coupling.to_string()));
        }
        write!(f, "Particle[idx:{} spin:{} Neighbours:({})]", self.id, self.spin, output)
    }
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn add_remove_neighbour(){
        let mut particles:Vec<Particle> = Vec::new();

        particles.push(Particle::new(
            0,
            Spin::Up));

        particles.push(Particle::new(
            1,
            Spin::Up));
                
        particles.push(Particle::new(
            2,
            Spin::Up));

        particles[0]
            .link_neighbour(Link::new(1,1.0))
            .link_neighbour(Link::new(2,1.0));
        
        let link:Link;
        {
            let tmp = particles[0].find_link(1).unwrap();
            link = tmp.clone()
        }
        particles[1]
            .link_neighbour(
                Link::new_from(
                    0,
                    &link));
                
        assert_eq!(particles[0].links[0].target_id, 1);

        particles[0]
            .remove_neighbour(1);

        assert_eq!(particles[0].links[0].target_id, 2);

        particles[0]
            .remove_neighbour(2)
            .remove_neighbour(9000);
        assert_eq!(particles[0].links.len(), 0);

    }
}
