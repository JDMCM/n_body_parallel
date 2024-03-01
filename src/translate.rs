
const PI: f64 = 3.141592653589793;
const SOLAR_MASS: f64 = 4.0 * PI * PI;
const YEAR: f64 = 365.24;
const N_BODIES: usize = 5;



#[derive(Clone, Copy)]
struct Nbody{
    x:f64,
    y:f64,
    z:f64,
    vx:f64,
    vy:f64,
    vz:f64,
    mass:f64,
}



pub fn circular_orbits(n: usize) -> Vec<Nbody> {
    let mut particle_buf = vec![];
      particle_buf.push(Nbody {
        x: 0.0,y: 0.0,z: 0.0,vx: 0.0,vy: 0.0,vz: 0.0,mass: 1.0*SOLAR_MASS,
      });
  
      for i in 0..n {
          let d = 0.1 + ((i as f64) * 5.0 / (n as f64));
          let v = f64::sqrt(1.0 / d);
          let theta = fastrand::f64() * 6.28;
          let x1 = d * f64::cos(theta);
          let y1 = d * f64::sin(theta);
          let vx1 = -v * f64::sin(theta);
          let vy1 = v * f64::cos(theta);
          particle_buf.push(Nbody {
              x: x1,y: y1,z: 0.0,vx: vx1,vy: vy1,vz: 0.0,
              mass: 1e-14*SOLAR_MASS,
             
          });
      }
      particle_buf
  }


fn step(&mut bodies: &mut Vec<Nbody>, dt: f64){
    for i in 0..bodies.len() {
        if i +1 < bodies.len() {
            for j in i+1..bodies.len(){
                // get the distance between the objects
                let dx: f64 = bodies[i].x - bodies[j].x;
                let dy: f64 = bodies[i].y - bodies[j].y;
                let dz: f64 = bodies[i].z - bodies[j].z;
                let d2: f64 = dx*dx + dy*dy + dz*dz;
                //get the magnitude of the force over a period of time 
                let magi:f64 = (dt*bodies[i].mass)/(d2*d2.sqrt());
                let magj:f64 = (dt*bodies[j].mass)/(d2*d2.sqrt());
                //update the velocities of the objects and proj mag onto eac distance comp
                bodies[i].vx += dx*magj;
                bodies[i].vy += dy*magj;
                bodies[i].vz += dz*magj;
                bodies[j].vx -= dz*magi;
                bodies[j].vy -= dz*magi;
                bodies[j].vz -= dz*magi;
            }
        }
    }
}

fn advance(bodies: &mut Vec<Nbody>, dt: f64, steps: i32) {
    for _ in 0..steps {
        step(bodies,dt)
    }
}

fn energy(bodies: &Vec<Nbody>) -> f64 {
    let mut e = 0.0;
    for i in 0..bodies.len() {
        if i +1 < bodies.len() {
            e += 0.5*bodies[i].mass*(bodies[i].vx*bodies[i].vx+bodies[i].vy*bodies[i].vy+bodies[i].vz*bodies[i].vz);
            for j in i+1..bodies.len(){
                // get the distance between the objects
                let dx: f64 = bodies[i].x - bodies[j].x;
                let dy: f64 = bodies[i].y - bodies[j].y;
                let dz: f64 = bodies[i].z - bodies[j].z;
                let d: f64 = (dx*dx + dy*dy + dz*dz).sqrt();
                e -= bodies[i].mass*bodies[j].mass / d;
            }
        }

    }
    return e;    
}

fn offset_momentum(bodies: &mut Vec<Nbody>){
    let mut px: f64 = 0.0;
    let mut py: f64 = 0.0;
    let mut pz: f64 = 0.0;
    for i in 1..bodies.len() {
        px += bodies[i].mass*bodies[i].vx;
        py += bodies[i].mass*bodies[i].vy;
        pz += bodies[i].mass*bodies[i].vz;
    }
    bodies[0].vx = px/bodies[0].mass;
    bodies[0].vy = py/bodies[0].mass;
    bodies[0].vz = pz/bodies[0].mass;
}

fn main() {
    let n = std::env::args_os().nth(1)
        .and_then(|s| s.into_string().ok())
        .and_then(|n| n.parse().ok())
        .unwrap_or(1000);
    let mut bodies = circular_orbits(10);

    offset_momentum(&mut bodies);
    println!("{:.9}", energy(&bodies));

    advance(&mut bodies, 0.01, n);

    println!("{:.9}", energy(&bodies));
}