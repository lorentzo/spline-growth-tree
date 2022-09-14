
#
# Blender 3.2.2/3.3.0
#

import bpy
import mathutils
import bmesh

import math

import numpy as np

# https://en.wikipedia.org/wiki/Helix
# Helix in x-y growing in z.
# t - parameter
# r - radius
# f - phase
# P - amount of squashness
def helix(t, r, P, f):
    x = r * math.cos(t + f)
    y = r * math.sin(t + f)
    z = P * t
    return mathutils.Vector((x,y,z))

def create_helix_points(r, P, f, length, resolution):
    points = []
    for t in np.linspace(0, length, resolution):
        point = helix(t, r, P, f)
        points.append(point)
    return points

def select_activate_only(objects=[]):
    for obj in bpy.data.objects:
        obj.select_set(False)
    bpy.context.view_layer.objects.active = None 
    for obj in objects:
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj

# https://blender.stackexchange.com/questions/220072/check-using-name-if-a-collection-exists-in-blend-is-linked-to-scene
def create_collection_if_not_exists(collection_name):
    if collection_name not in bpy.data.collections:
        new_collection = bpy.data.collections.new(collection_name)
        bpy.context.scene.collection.children.link(new_collection) #Creates a new collection

def create_icosphere(radius=1.0, name="ico", collection_name=None):
    bm = bmesh.new()
    # Create icosphere.
    # https://docs.blender.org/api/current/bmesh.ops.html#bmesh.ops.create_icosphere
    bmesh.ops.create_icosphere(bm, subdivisions=1, radius=radius, matrix=mathutils.Matrix.Identity(4), calc_uvs=False)
    object_mesh = bpy.data.meshes.new(name+"_mesh")
    bm.to_mesh(object_mesh)
    obj = bpy.data.objects.new(name+"_obj", object_mesh)
    if collection_name == None:
        bpy.context.collection.objects.link(obj)
    else:
        create_collection_if_not_exists(collection_name)
        bpy.data.collections[collection_name].objects.link(obj)
    bm.free()
    return obj

def create_point_light(name="plight", color=(1, 0.376144, 0.0154181), energy=50, collection_name=None):
    point_light_data = bpy.data.lights.new(name=name+"_data", type="POINT")
    point_light_object = bpy.data.objects.new(name=name+"_object", object_data=point_light_data)
    point_light_object.data.color = color
    point_light_object.data.energy = energy
    if collection_name == None:
        bpy.context.collection.objects.link(point_light_object)
    else:
        create_collection_if_not_exists(collection_name)
        bpy.data.collections[collection_name].objects.link(point_light_object)
    return point_light_object


def create_instance(base_obj,
                    translate=mathutils.Vector((0,0,0)), 
                    scale=1.0,
                    rotate=("Z", 0.0),
                    basis=mathutils.Matrix.Identity(4),
                    tbn=mathutils.Matrix.Identity(4),
                    collection_name=None):
    # Create instance.
    inst_obj = bpy.data.objects.new(name=base_obj.name+"_inst", object_data=base_obj.data)
    #inst_obj = base_obj.copy()
    #inst_obj.data = base_obj.data.copy()
    #inst_obj.animation_data_clear()
    #inst_obj.location = mathutils.Vector((0,0,0))
    # Perform translation, rotation, scaling and moving to target coord system for instance.
    mat_rot = mathutils.Matrix.Rotation(math.radians(rotate[1]), 4, rotate[0])
    mat_trans = mathutils.Matrix.Translation(translate)
    mat_sca = mathutils.Matrix.Scale(scale, 4) # TODO: figure out how to scale in given vector direction
    # TODO: If I am using `tbn` as basis then it sould go last, If I use `matrix_basis` as basis then it should go first.
    # `tbn` matrix is usually constructed for samples on base geometry using triangle normal. Therefore, it only contains
    # information about rotation.
    inst_obj.matrix_basis = basis @ mat_trans @ mat_rot @ mat_sca @ tbn  # TODO: is matrix_basis correct to be used for this?
    # Store to collection.
    if collection_name == None:
        bpy.context.collection.objects.link(inst_obj)
    else:
        create_collection_if_not_exists(collection_name)
        bpy.data.collections[collection_name].objects.link(inst_obj)
    return inst_obj

def point_on_object(obj, point):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bvh = mathutils.bvhtree.BVHTree.FromBMesh(bm)
    nearest = bvh.find_nearest(point)
    bm.free()
    return nearest[0]

def project_points_on_object(base_obj=None, points=[]):
    projected_points = []
    bm = bmesh.new()
    bm.from_mesh(base_obj.data)
    bvh = mathutils.bvhtree.BVHTree.FromBMesh(bm)
    for point in points:
        nearest = bvh.find_nearest(point)
        projected_points.append(mathutils.Vector(nearest[0])) # take only position
    bm.free()
    return projected_points


def interpolate_two_points(starting_point=mathutils.Vector((0,0,0)), ending_point=mathutils.Vector((0,0,0)), n_subdivisions=0):
    points_dists = []
    t = 0
    dt = 1.0
    if n_subdivisions > 0:
        dt = 1.0 / (n_subdivisions+1)
    for subdiv_i in range(n_subdivisions+2):
        pi = (1.0 - t) * ending_point + t * starting_point
        points_dists.append((mathutils.Vector(pi), mathutils.Vector(starting_point-pi).length))
        t = t + dt
    points_dists.sort(key=lambda x:x[1])
    points = []
    for point_dist in points_dists:
        points.append(point_dist[0])
    return points

def subsample_points(points=[], min_dist=100.0):
    min_elem_dist = min_dist
    subsampled_points = []
    for i in range(len(points)-1):
        p1 = mathutils.Vector(points[i])
        p2 = mathutils.Vector(points[i+1])
        p1_p2_vec = mathutils.Vector(p1 - p2)
        p1_p2_len = p1_p2_vec.length
        n_subdiv = int(math.floor(p1_p2_len / min_elem_dist))
        if n_subdiv > 0:
            interpolated_points = interpolate_two_points(starting_point=p1, ending_point=p2, n_subdivisions=n_subdiv)
            for interpolated_point in interpolated_points:
                subsampled_points.append(interpolated_point)
        else:
            subsampled_points.append(p1)
            subsampled_points.append(p2)
    return subsampled_points


# https://blenderscripting.blogspot.com/2011/05/blender-25-python-bezier-from-list-of.html
def create_curve(name="curve", points=[], bevel=0, bevel_start=0.0, bevel_end=0.0, collection_name=None):
    # Create curve data.
    curve_data = bpy.data.curves.new(name=name, type="CURVE")
    curve_data.dimensions = "3D"
    curve_data.bevel_depth = bevel
    curve_data.bevel_factor_start = bevel_start
    curve_data.bevel_factor_end = bevel_end
    # Populate curve data.
    poly_curve = curve_data.splines.new("POLY")
    poly_curve.points.add(len(points)-1)
    for pi in range(len(points)):
        poly_curve.points[pi].co = (points[pi][0], points[pi][1], points[pi][2], 1)
    # Create curve object from curve data and add it to the scene.
    curve_object = bpy.data.objects.new(name=name+"_obj", object_data=curve_data)
    if collection_name == None:
        bpy.context.collection.objects.link(curve_object)
    else:
        create_collection_if_not_exists(collection_name)
        bpy.data.collections[collection_name].objects.link(curve_object)
    return curve_object

# Option1: https://docs.blender.org/api/current/bpy.ops.curve.html#bpy.ops.curve.smooth
def smooth_curve1(curve_obj=None, smooth_iterations=5):
    select_activate_only(objects=[curve_obj])
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.curve.select_all(action='SELECT')
    for i in range(smooth_iterations):
        bpy.ops.curve.smooth()
    bpy.ops.object.mode_set(mode='OBJECT')



# Option2: https://docs.blender.org/api/current/bpy.types.SmoothModifier.html
def smooth_curve2(curve_obj=None, smooth_factor=0.5, smooth_iterations=5):
    smooth_modifier = curve_obj.modifiers.new(name="SMOOTH", type="SMOOTH")
    smooth_modifier.factor = smooth_factor
    smooth_modifier.iterations = smooth_iterations
    select_activate_only(objects=[curve_obj])
    #bpy.ops.object.modifier_apply(modifier="SMOOTH")

def jitter_points(points=[], jitter_amount=1.0, jitter_freq=1.0):
    jittered_points = []
    for point in points:
        jittered_point = mathutils.Vector(point + mathutils.noise.noise_vector(point * jitter_freq) * jitter_amount)
        jittered_points.append(jittered_point)
    return jittered_points

def curl_points(points=[], curl_scale=1.0, freq_scale=1.0):
    curled_points = []
    for point in points:
        octaves = 4
        hard = False
        curled_point = point + mathutils.noise.turbulence_vector(point, octaves, hard, noise_basis='BLENDER', amplitude_scale=curl_scale, frequency_scale=freq_scale)
        curled_points.append(curled_point)
    return curled_points

# https://vividfax.github.io/2021/01/14/blender-materials.html
def create_and_assign_emissive_material(base_obj=None, mat_name="emissive", color=(1,1,1,1), energy=1.0):
    base_obj.data.materials.clear()
    emission = bpy.data.materials.new(name=mat_name)
    emission.use_nodes = True
    emission_nodes = emission.node_tree.nodes
    emission_links = emission.node_tree.links
    emission_nodes.clear()
    emission_links.clear()
    emission_output = emission_nodes.new(type='ShaderNodeOutputMaterial')
    emission_shader = emission_nodes.new(type="ShaderNodeEmission")
    emission_shader.inputs[0].default_value = color
    emission_shader.inputs[1].default_value = energy
    emission_links.new(emission_shader.outputs[0], emission_output.inputs[0])
    base_obj.data.materials.append(emission)



def main():
    # Get the base object.
    tree = bpy.data.objects["base_obj"]

    # get the spline desribing the animation path.
    # NOTE for creating animation_paths: create curve path (nurbs path), delete all points except one. Extrude and transform from that one.
    # This is needed since the animation direction will be correct if this is done.
    for animation_path_object in bpy.data.collections['animation_paths'].objects:

        animation_path_spline = animation_path_object.data.splines[0]
        animation_path_points = animation_path_spline.points
        animation_path_point_coords = []
        for animation_path_point in animation_path_points:
            coord = mathutils.Vector((animation_path_point.co[0], animation_path_point.co[1], animation_path_point.co[2]))
            animation_path_point_coords.append(coord)

        # Project animation path points on base object.
        animation_path_points_projected_on_base_object = project_points_on_object(base_obj=tree, points=animation_path_point_coords)
        # Subsample projected points if needed.
        animation_path_points_projected_on_base_object_subsampled = subsample_points(points=animation_path_points_projected_on_base_object, min_dist=1.0)
        
        # Create multiple animation paths by jittering input points.
        n_animation_paths = 10
        for i in range(n_animation_paths):

            # Jitter points.
            jitter_amount_max = 2.0
            jitter_amount_min = 1.0
            jitter_amount_t = mathutils.noise.random()
            jitter_amount = (1.0 - jitter_amount_t) * jitter_amount_min + jitter_amount_t * jitter_amount_max
            jitter_freq_max = 3.0
            jitter_freq_min = 1.0
            jitter_freq_t = mathutils.noise.random()
            jitter_freq = (1.0 - jitter_freq_t) * jitter_freq_min + jitter_freq_t * jitter_freq_max
            animation_path_points_projected_on_base_object_subsampled_jittered = jitter_points(points=animation_path_points_projected_on_base_object_subsampled, jitter_amount=jitter_amount, jitter_freq=jitter_freq)
            
            # Randomize number of animation path points.
            animation_path_points_max = len(animation_path_points_projected_on_base_object_subsampled_jittered)
            animation_path_points_min = len(animation_path_points_projected_on_base_object_subsampled_jittered) / 4.0
            animation_path_points_t = mathutils.noise.random()
            animation_path_points_n = int(math.ceil((1.0 - animation_path_points_t) * animation_path_points_min + animation_path_points_t * animation_path_points_max))
            for i in range(animation_path_points_max - animation_path_points_n):
                animation_path_points_projected_on_base_object_subsampled_jittered.pop()
            
            # Add curling top.
            # Take last two points to create a vector, create a new point in direction of that vector, subsample and curl.
            ultimate_point = animation_path_points_projected_on_base_object_subsampled_jittered[-1]
            penultimate_point = animation_path_points_projected_on_base_object_subsampled_jittered[-2]
            growth_dir = mathutils.Vector(ultimate_point - penultimate_point).normalized()
            print(growth_dir)
            curling_top_point_offset = 15.0
            curling_top_point = ultimate_point + growth_dir * curling_top_point_offset
            curling_points = subsample_points(points=[ultimate_point, curling_top_point], min_dist=1.0)
            curling_points = curl_points(points=curling_points, curl_scale=1.0, freq_scale=1.0)
            animation_path_points_projected_on_base_object_subsampled_jittered.extend(curling_points)

            # Create curve using generated points.
            growth_curve = create_curve(name="tree_spline", points=animation_path_points_projected_on_base_object_subsampled_jittered, bevel=0.0, bevel_start=0.0, bevel_end=0.0, collection_name="procedural_growth")
            smooth_curve1(curve_obj=growth_curve, smooth_iterations=15)
            
            # Create taper curve to control curve shape.
            taper_curve_points = [mathutils.Vector((-2,0,0)), mathutils.Vector((-1,0,0)), mathutils.Vector((0,0,0)), mathutils.Vector((1,0,0)), mathutils.Vector((2,0,0))]
            taper_curve = create_curve(name="taper_curve", points=taper_curve_points, bevel=0.0, bevel_start=0.0, bevel_end=0.0, collection_name="procedural_growth")
            growth_curve.data.taper_object = taper_curve
            for p in taper_curve.data.splines[0].points:
                p.co[1] = 0.1
            taper_curve.data.splines[0].points[0].co[1] += 0.4
            taper_curve.data.splines[0].points[1].co[1] += 0.15
            taper_curve.data.splines[0].points[2].co[1] += 0.09
            taper_curve.data.splines[0].points[3].co[1] += 0.05
            
            # Animate growth.
            delta_keypoint_frame = 20
            n_growth_animation_frames = 300
            bevel_depth_delta = 0.4 / n_growth_animation_frames
            bevel_factor_delta = 1.0 / n_growth_animation_frames
            for frame in range(n_growth_animation_frames+1):
                if frame % delta_keypoint_frame == 0:
                    growth_curve.data.keyframe_insert(data_path="bevel_depth", frame=frame)
                    growth_curve.data.keyframe_insert(data_path="bevel_factor_start", frame=frame)
                growth_curve.data.bevel_depth += bevel_depth_delta
                growth_curve.data.bevel_factor_start += bevel_factor_delta

            # Animate movement of splines.
            # TODO.

            # Spawn fliers moving alongside animation path.
            # TODO.

            if mathutils.noise.random() > 0.5:
                # Animate lights inside curves.
                light = create_point_light(name="plight", color=(1, 0.376144, 0.0154181), energy=50, collection_name="lights")
                n_light_steps = 200
                starting_light_animation_frame = n_growth_animation_frames
                delta_keypoint_light_frame = int(math.ceil(n_light_steps / len(animation_path_points_projected_on_base_object_subsampled_jittered)))
                light_keyframe = starting_light_animation_frame
                for point in growth_curve.data.splines[0].points:
                    p = mathutils.Vector((point.co[0], point.co[1], point.co[2]))
                    light.location = p
                    light.keyframe_insert(data_path="location", frame=light_keyframe)
                    light_keyframe += delta_keypoint_light_frame
                # Animate spheres flying in the sky and staying there.
                starting_stars_animation_frame = starting_light_animation_frame + n_light_steps
                p1 = growth_curve.data.splines[0].points[-2].co
                p1_loc = mathutils.Vector((p1[0], p1[1], p1[2]))
                p2 = growth_curve.data.splines[0].points[-1].co
                p2_loc = mathutils.Vector((p2[0], p2[1], p2[2]))
                p1p2_dir = mathutils.Vector(p2_loc - p1_loc).normalized()
                star = create_icosphere(radius=1.0, name="star", collection_name="stars")
                create_and_assign_emissive_material(base_obj=star, mat_name="emissive", color=(1, 0.385386, 0.0399476, 1), energy=50.0)
                star.scale = mathutils.Vector((0.0, 0.0, 0.0))
                star.keyframe_insert(data_path="scale", frame=0)
                star.keyframe_insert(data_path="scale", frame=starting_stars_animation_frame-1)
                star.location = p2_loc
                star.scale = mathutils.Vector((0.3, 0.3, 0.3))
                star.keyframe_insert(data_path="location", frame=starting_stars_animation_frame)
                star.keyframe_insert(data_path="scale", frame=starting_stars_animation_frame)
                n_star_steps = 30
                stars_min_dist = 150.0
                stars_max_dist = 20.0
                stars_dist_t = mathutils.noise.random()
                star_dist_curr = (1 - stars_dist_t) * stars_min_dist + stars_dist_t * stars_max_dist
                star_movement_delta = star_dist_curr / n_star_steps
                star_location_start = p2_loc
                for frame in range(n_star_steps):
                    star.keyframe_insert(data_path="location", frame=starting_stars_animation_frame + frame)
                    star.location = star.location + p1p2_dir * star_movement_delta
                # TODO: noise movement
                # TODO: noise emission
                # TODO: add emissive material
                    
#
# Script entry point.
#
if __name__ == "__main__":
    main()