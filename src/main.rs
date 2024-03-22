use rayon::prelude::*;

const PI: f64 = 3.141592653589793;
const SOLAR_MASS: f64 = 4.0 * PI * PI;
const YEAR: f64 = 365.24;
const N_BODIES: usize = 5;



#[derive(Clone, Copy, Debug)]
struct Nbody{
    x:f64,
    y:f64,
    z:f64,
    vx:f64,
    vy:f64,
    vz:f64,
    mass:f64,
    
}





fn create_testbodies() ->Vec<Nbody> {
let mut testbodies:Vec<Nbody> = vec![];

testbodies.push(Nbody {
    x: 0.0, y: 0.0, z: 0.0,
    vx: 0.0, vy: 0.0, vz: 0.0,
    mass: SOLAR_MASS, 
});
// Jupiter
testbodies.push(Nbody {
    x: 4.84143144246472090e+00,
    y: -1.16032004402742839e+00,
    z: -1.03622044471123109e-01,
    vx: 1.66007664274403694e-03 * YEAR,
    vy: 7.69901118419740425e-03 * YEAR,
    vz: -6.90460016972063023e-05 * YEAR,
    mass: 9.54791938424326609e-04 * SOLAR_MASS,
    
});
// Saturn
testbodies.push(Nbody {
    x: 8.34336671824457987e+00,
    y: 4.12479856412430479e+00,
    z: -4.03523417114321381e-01,
    vx: -2.76742510726862411e-03 * YEAR,
    vy: 4.99852801234917238e-03 * YEAR,
    vz: 2.30417297573763929e-05 * YEAR,
    mass: 2.85885980666130812e-04 * SOLAR_MASS,
    
});
// Uranus
testbodies.push(Nbody {
    x: 1.28943695621391310e+01,
    y: -1.51111514016986312e+01,
    z: -2.23307578892655734e-01,
    vx: 2.96460137564761618e-03 * YEAR,
    vy: 2.37847173959480950e-03 * YEAR,
    vz: -2.96589568540237556e-05 * YEAR,
    mass: 4.36624404335156298e-05 * SOLAR_MASS,
    
});
// Neptune
testbodies.push(Nbody {
    x: 1.53796971148509165e+01,
    y: -2.59193146099879641e+01,
    z: 1.79258772950371181e-01,
    vx: 2.68067772490389322e-03 * YEAR,
    vy: 1.62824170038242295e-03 * YEAR,
    vz: -9.51592254519715870e-05 * YEAR,
    mass: 5.15138902046611451e-05 * SOLAR_MASS,
    
});

return testbodies
}


fn step(bodies: &mut Vec<Nbody>, dt: f64){
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
                bodies[i].vx -= dx*magj;
                bodies[i].vy -= dy*magj;
                bodies[i].vz -= dz*magj;
                bodies[j].vx += dx*magi;
                bodies[j].vy += dy*magi;
                bodies[j].vz += dz*magi;
            }
        }


    }
    for i in 0..bodies.len() {
        bodies[i].x += bodies[i].vx*dt;
        bodies[i].y += bodies[i].vy*dt;
        bodies[i].z += bodies[i].vz*dt;
    }
}

fn advance(bodies: &mut Vec<Nbody>, dt: f64, steps: i32) {
    for _ in 0..steps {
        step(bodies,dt)
    }
}

fn make_accels(bodies: & Vec<Nbody>,dt: f64)->Vec<(f64,f64,f64)>{
    let new_vels = bodies.into_par_iter().map(|i| {
        let mut new_vel: (f64,f64,f64) = (i.vx,i.vy,i.vz);
        
        for j in  bodies {
            let dx: f64 = i.x - j.x;
            let dy: f64 = i.y - j.y;
            let dz: f64 = i.z - j.z;
            let d2: f64 = dx*dx + dy*dy + dz*dz;
            if d2 != 0.0 {
                let magi:f64 = (dt*j.mass)/(d2*d2.sqrt());
                new_vel.0 -= dx*magi;
                new_vel.1 -= dy*magi;
                new_vel.2 -= dz*magi;

            }
        }
        return new_vel;
    }).collect();
    return new_vels;
}

fn pstep(bodies: &mut Vec<Nbody>,dt: f64) {
    let mut new_vel = make_accels(bodies, dt);
    bodies.par_iter_mut().zip(&mut new_vel).for_each(|(i,j)| {
        i.vx = j.0;
        i.vy = j.1;
        i.vz = j.2;
        i.x  += j.0*dt;
        i.y  += j.1*dt;
        i.z  += j.2*dt;

    })

    
    // for i in 0..bodies.len() {
    //     bodies[i].vx = new_vel[i].0;
    //     bodies[i].vy = new_vel[i].1;
    //     bodies[i].vz = new_vel[i].2;
    //     bodies[i].x +=  bodies[i].vx*dt;
    //     bodies[i].y +=  bodies[i].vy*dt;
    //     bodies[i].z +=  bodies[i].vz*dt;
    // }

}

fn padvance(bodies: &mut Vec<Nbody>, dt: f64, steps: i32){
    for _ in 0..steps {
        pstep(bodies,dt);
        //println!("{:?},{:?},{:?}", bodies[1].x ,bodies[1].y, bodies[1].z);
    }
}


    


fn worse_step(bd: &mut Vec<Nbody>, dt: f64) {
    let mut mtopx:f64 = 0.0;
    let mut mtopy:f64 = 0.0;
    let mut mtopz:f64 = 0.0;
    let mut mbot:f64  = 0.0;
    //get the center mass of the total system
    for i in bd.iter() {
        mtopx += i.x*i.mass;
        mtopy += i.y*i.mass;
        mtopz += i.z*i.mass;
        mbot  += i.mass;
    }
    //get the center of mass of the entire system except
    //the part we are finding the accelleration on
    for i in bd.iter_mut() {
       let m_sys = mbot -i.mass;
       let r_x = i.x-((mtopx - i.x*i.mass)/m_sys);
       let r_y = i.y-((mtopy - i.y*i.mass)/m_sys);
       let r_z = i.z-((mtopz - i.z*i.mass)/m_sys);
       let r2 = r_x*r_x+r_y*r_y+r_z*r_z;
       let mag_a = (dt*m_sys)/(r2*r2.sqrt());
       i.vx -= r_x*mag_a;
       i.vx -= r_y*mag_a;
       i.vx -= r_z*mag_a;
       i.x += i.vx*dt;
       i.y += i.vy*dt;
       i.z += i.vz*dt;

    }
}

fn p_energy(bodies: & Vec<Nbody>)->f64{
    let e_s = (0..bodies.len()).into_par_iter().map(|i| {
        let mut e: f64 = 0.5*bodies[i].mass*(bodies[i].vx*bodies[i].vx+bodies[i].vy*bodies[i].vy+bodies[i].vz*bodies[i].vz);
        for j in  i+1..bodies.len() {
             // get the distance between the objects
             let dx: f64 = bodies[i].x - bodies[j].x;
             let dy: f64 = bodies[i].y - bodies[j].y;
             let dz: f64 = bodies[i].z - bodies[j].z;
             let d: f64 = (dx*dx + dy*dy + dz*dz).sqrt();
             //println!("{} {} {} {}",i,j,d,e);
             e -= bodies[i].mass*bodies[j].mass / d;
        }
        
        return e;
    });
    let pot: f64 = e_s.sum();
    
    return pot;
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
                //println!("{} {} {} {}",i,j,d,e);
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
    bodies[0].vx = -px/bodies[0].mass;
    bodies[0].vy = -py/bodies[0].mass;
    bodies[0].vz = -pz/bodies[0].mass;
}

pub fn circular_orbits(n: usize) -> Vec<Nbody> {
    let mut particle_buf = vec![];
      particle_buf.push(Nbody {
        x: 0.0,y: 0.0,z: 0.0,vx: 0.0,vy: 0.0,vz: 0.0,mass: 1.0*SOLAR_MASS, 
      });
  
      for i in 1..n+1 {
          let d = 0.1 + ((i as f64) * 5.0 / (n as f64));
          let v = f64::sqrt(1.0 / d);
          let theta = fastrand::f64() * 6.28;
          let x1 = d * f64::cos(theta);
          let y1 = d * f64::sin(theta);
          let vx1 = -v * f64::sin(theta)*365.0;
          let vy1 = v * f64::cos(theta)*365.0;
          particle_buf.push(Nbody {
              x: x1,y: y1,z: 0.0,vx: vx1,vy: vy1,vz: 0.0,
              mass: 1e-14*SOLAR_MASS, 
             
          });
      }
      particle_buf
}

fn main() {
    let n = std::env::args_os().nth(1)
        .and_then(|s| s.into_string().ok())
        .and_then(|n| n.parse().ok())
        .unwrap_or(1000);
    let mut bodies = circular_orbits(1000);
    let mut bodies1 = circular_orbits(1000);
    //create_testbodies(); //circular_orbits(10);
    //println!("{:?}", bodies);


    offset_momentum(&mut bodies);
    offset_momentum(&mut bodies1);
    //parallel
    println!("{}", p_energy(&bodies));


    padvance(&mut bodies, 0.01,n);

    println!("{}", p_energy(&bodies));
    //seq
    println!("{}", energy(&bodies1));


    advance(&mut bodies1, 0.01,n);

    println!("{}", energy(&bodies1));
}