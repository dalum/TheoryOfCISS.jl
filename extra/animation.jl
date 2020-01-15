euler_angles(p) = (n = norm(p); [0.0, acos(-p[3] / n), atan(-p[2], -p[1])])

# Find the actual crossing point
x0, y0 = TheoryOfCISS.find_crossing(
    1.23,
    2.39,
    14,
    threshold=1e-6,
    attenuation=1e-2
);

# Calculate the D and C vectors
vs = TheoryOfCISS.atpoint(
    x -> x[:v],
    x0,
    y0,
    N = 14,
    calc = TheoryOfCISS.calc_data1
);

D = 2real.(vs)
C = 2imag.(vs)

# Calculate the contribution to Λ12 from each atom
vxs, vys, vzs = [
    TheoryOfCISS.atpoint(
        x -> (xx -> xx[i]).(x[:vs]),
        x0 - 0.05,
        y0,
        N = 14,
        calc = TheoryOfCISS.calc_data2
    )
    for i in 1:3
];

ds = dxs, dys, dzs = real.(vxs), real.(vys), real.(vzs)
cs = cxs, cys, czs = imag.(vxs), imag.(vys), imag.(vzs)

λ = 8e-3

# Normalize using the total sum of contributions
N_d = sqrt(sum(dxs)^2 + sum(dys)^2 + sum(dzs)^2);
N_c = sqrt(sum(cxs)^2 + sum(cys)^2 + sum(czs)^2);
# Cumulatively add up contributions
cumulated_d = collect(zip(cumsum(dxs), cumsum(dys), cumsum(dzs)))
cumulated_c = collect(zip(cumsum(cxs), cumsum(cys), cumsum(czs)))

# Compute the scaling and Euler angles
scales = norm.(polarizations)
angles = euler_angles.(polarizations)

script = """
import bpy

electron_vector = bpy.context.scene.objects["ElectronVector"]

scales = [$(join(scales, ","))]
angles = [$(join(angles, ","))]

bpy.ops.screen.frame_jump(end=False)

frame = 10

for (scale, (x, y, z)) in zip(scales, angles):
    bpy.context.scene.frame_set(frame)
    electron_vector.scale[2] = scale
    electron_vector.rotation_euler[0] = x
    electron_vector.rotation_euler[1] = y
    electron_vector.rotation_euler[2] = z
    bpy.ops.anim.keyframe_insert()
    frame += 10

bpy.ops.screen.frame_jump(end=False)
"""

write("animate.py", script)
