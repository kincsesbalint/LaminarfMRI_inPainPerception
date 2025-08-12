import nibabel as nib
import numpy as np
import trimesh

def load_freesurfer_surface(surface_file):
    """Load a FreeSurfer surface file (.pial or .white) using nibabel."""
    surf = nib.freesurfer.read_geometry(surface_file)
    vertices, faces = surf[0], surf[1]
    return vertices, faces

def compute_vertex_normal(vertices, faces, vertex_index):
    """Compute the normal at a given vertex using neighboring face normals."""
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    normals = mesh.vertex_normals
    return normals[vertex_index]

def find_pial_intersection(wm_vertex, normal, pial_vertices, pial_faces):
    """Find intersection of the normal ray with the pial surface."""
    mesh = trimesh.Trimesh(vertices=pial_vertices, faces=pial_faces, process=False)
    ray_origins = np.array([wm_vertex])
    ray_directions = np.array([normal])

    locations, index_ray, index_tri = mesh.ray.intersects_location(ray_origins, ray_directions)
    
    if len(locations) > 0:
        # Take the closest intersection point
        distances = np.linalg.norm(locations - wm_vertex, axis=1)
        intersection = locations[np.argmin(distances)]
        return intersection
    else:
        print("No intersection found with the pial surface.")
        return None

def sample_along_line(start, end, num_samples=10):
    """Generate points along the line between start and end."""
    return np.linspace(start, end, num_samples)

# Paths to FreeSurfer surface files
wm_surface_file = "subject/surf/lh.white"  # White matter surface
pial_surface_file = "subject/surf/lh.pial"  # Pial surface

# Load surfaces
wm_vertices, wm_faces = load_freesurfer_surface(wm_surface_file)
pial_vertices, pial_faces = load_freesurfer_surface(pial_surface_file)

# Choose a white matter vertex index (e.g., 10000)
vertex_index = 10000
wm_vertex = wm_vertices[vertex_index]

# Compute normal at the white matter vertex
normal = compute_vertex_normal(wm_vertices, wm_faces, vertex_index)

# Find intersection with the pial surface along the normal direction
pial_intersection = find_pial_intersection(wm_vertex, normal, pial_vertices, pial_faces)

if pial_intersection is not None:
    # Sample along the perpendicular line
    sampled_points = sample_along_line(wm_vertex, pial_intersection, num_samples=10)

    print("Sampled Points Along the Perpendicular Line:")
    print(sampled_points)
