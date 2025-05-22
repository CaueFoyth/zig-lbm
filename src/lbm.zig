const VelSet: type = enum{
    D2Q9,
};

const vel_set_use: VelSet = VelSet.D2Q9;

//Dimensões
const dim = switch(vel_set_use){
    .D2Q9 => 2,
};

//Número de populações
const n_pop = switch(vel_set_use){
    .D2Q9 => 9,
};

//Direções das populações
const pop_dir: [n_pop][dim]u8 = switch(vel_set_use){
    .D2Q9 => [_][dim]u8{[_]u8{0, 0}, [_]u8{1, 0}, [_]u8{0, 1}, [_]u8{-1, 0}, [_]u8{0, -1}, [_]u8{1, 1}, [_]u8{-1, 1}, [_]u8{-1, -1}, [_]u8{1, -1}},
};

//Pesos das populações
const pop_weights: [n_pop]f32 = switch(vel_set_use){
    .D2Q9 => [_]f32{4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0},
};

const domain_size: [dim]u32 = .{32,32};

fn idx2pos(idx: usize) [dim]u32 {
    if (dim == 2){
        return .{idx % domain_size[0], idx / domain_size[0]};
    } else{
        return .{idx % domain_size[0], (idx / domain_size[0]) % domain_size[1], idx % (domain_size[0] * domain_size[1])};
    }
}

fn pos2idx(pos: [dim]u32) usize {
    if (dim == 2){
        return pos[0] + pos[1] * domain_size[0];
    }else{
        return pos[0] + domain_size[0] * (pos[1] + pos[2] * domain_size[1]);
    }
}

const cs2: f32 = 1.0 / 3.0;
const tau: f32 = 0.9;

fn dot_prod(x: []f32, y: []f32) f32 {
    var sum: f32 = 0;
    for(x, y) |i, j| {
        sum += i * j;
    }
    return sum;
}
fn func_feq(rho: f32, u: [dim]f32, comptime i: usize) f32 {
    const uc = dot_prod(u, pop_dir[i]);
    const uu = dot_prod(u, pop_dir[i]);

    return rho * pop_weights[i] * (
        1 + uc / cs2 + uc * uc / (2 * cs2 * cs2) - uu / (2 * cs2)
    );
}

pub fn macroscopics(
    pop_arr: [][n_pop]f32, 
    rho_arr: []f32,
    u_arr: [][dim]f32
) void{
    var i: usize = 0;
    while(i < pop_arr.len) : (i += 1){
        const pop = pop_arr[i];
        var rho: f32 = 0;
        for (pop) |p| {
            rho += p;
        }
        var u: [dim]f32 = .{0} ** dim;
        for (pop, 0..) |p, j| {
            for (dim) |d|{
                u[d] += p * pop_dir[j];
            }
        }

        rho_arr[i] = rho;
        u_arr[i] = u;
    }
}

pub fn collision(
    pop_arr: [][n_pop]f32,
    rho_arr: []f32,
    u_arr: [][dim]f32
) void{
    for (pop_arr, 0..) |pop, idx| {
        const rho = rho_arr[idx];
        const u = u_arr[idx];
        for (pop, 0..) |f, i| {
            const feq = func_feq(rho, u, i);
            const f_col = (f - feq) / tau;
            pop[idx][i] = f_col;
        }
    }
}

pub fn streaming(popA_arr: [][n_pop]f32, popB_arr: [][n_pop]f32) void{
    for(popA_arr, 0..) |pop, idx| {
        const pos = idx2pos(idx);
        for (pop, 0..) |p, i|{
            const postTo: [dim]i32 = .{(@as(i32, pos[0]) + @as(i32, pop_dir[0]) + domain_size[0]) % domain_size[0], (@as(i32, pos[1]) + @as(i32, pop_dir[1]) + domain_size[1]) % domain_size[1]};
            const idxTo: usize = pos2idx(@as([2]32, postTo));
            popB_arr[idxTo][i] = p;
        } 
    }
}

pub fn run_time_step(
    popA_arr: [][n_pop]f32,
    popB_arr: [][n_pop]f32,
    rho_arr: []f32,
    u_arr: [][dim]f32,
    time_step: u32
) void{
    const popMain_arr = if (time_step % 2 == 0) popA_arr else popB_arr;
    const popAux_arr = if (time_step % 2 == 1) popA_arr else popB_arr;

    macroscopics(popMain_arr, rho_arr, u_arr);
    collision(popMain_arr, rho_arr, u_arr);
    streaming(popMain_arr, popAux_arr);
}