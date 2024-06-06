#This header file imports open3d with the correct version to support version 0.5.0 and 0.9.0
import open3d as o3d
from packaging import version
if version.parse('0.12.0') <= version.parse(o3d.__version__):
    PointCloud = o3d.geometry.PointCloud
    Vector3dVector = o3d.utility.Vector3dVector
    Feature = o3d.pipelines.registration.Feature
    read_point_cloud = o3d.io.read_point_cloud
    read_triangle_mesh = o3d.io.read_triangle_mesh
    registration_ransac_based_on_feature_matching = o3d.pipelines.registration.registration_ransac_based_on_feature_matching
    registration_icp = o3d.pipelines.registration.registration_icp
    TransformationEstimationPointToPoint = o3d.pipelines.registration.TransformationEstimationPointToPoint
    CorrespondenceCheckerBasedOnEdgeLength = o3d.pipelines.registration.CorrespondenceCheckerBasedOnEdgeLength
    CorrespondenceCheckerBasedOnDistance = o3d.pipelines.registration.CorrespondenceCheckerBasedOnDistance
    CorrespondenceCheckerBasedOnNormal = o3d.pipelines.registration.CorrespondenceCheckerBasedOnNormal
    TransformationEstimationPointToPlane = o3d.pipelines.registration.TransformationEstimationPointToPlane
    RANSACConvergenceCriteria = o3d.pipelines.registration.RANSACConvergenceCriteria
    KDTreeFlann = o3d.geometry.KDTreeFlann
elif version.parse('0.6.0') < version.parse(o3d.__version__):
    PointCloud = o3d.geometry.PointCloud
    Vector3dVector = o3d.utility.Vector3dVector
    Feature = o3d.registration.Feature
    read_point_cloud = o3d.io.read_point_cloud
    read_triangle_mesh = o3d.io.read_triangle_mesh
    registration_ransac_based_on_feature_matching = o3d.registration.registration_ransac_based_on_feature_matching
    registration_icp = o3d.registration.registration_icp
    TransformationEstimationPointToPoint = o3d.registration.TransformationEstimationPointToPoint
    CorrespondenceCheckerBasedOnEdgeLength = o3d.registration.CorrespondenceCheckerBasedOnEdgeLength
    CorrespondenceCheckerBasedOnDistance = o3d.registration.CorrespondenceCheckerBasedOnDistance
    CorrespondenceCheckerBasedOnNormal = o3d.registration.CorrespondenceCheckerBasedOnNormal
    TransformationEstimationPointToPlane = o3d.registration.TransformationEstimationPointToPlane
    RANSACConvergenceCriteria = o3d.registration.RANSACConvergenceCriteria
    KDTreeFlann = o3d.geometry.KDTreeFlann
else:
    PointCloud = o3d.PointCloud
    Vector3dVector = o3d.Vector3dVector
    Feature = o3d.Feature
    read_triangle_mesh = o3d.read_triangle_mesh
    read_point_cloud = o3d.read_point_cloud
    registration_ransac_based_on_feature_matching = o3d.registration_ransac_based_on_feature_matching
    TransformationEstimationPointToPoint = o3d.TransformationEstimationPointToPoint
    CorrespondenceCheckerBasedOnEdgeLength = o3d.CorrespondenceCheckerBasedOnEdgeLength
    registration_icp = o3d.registration_icp
    CorrespondenceCheckerBasedOnDistance = o3d.CorrespondenceCheckerBasedOnDistance
    CorrespondenceCheckerBasedOnNormal = o3d.CorrespondenceCheckerBasedOnNormal
    TransformationEstimationPointToPlane = o3d.TransformationEstimationPointToPlane
    RANSACConvergenceCriteria = o3d.RANSACConvergenceCriteria
    KDTreeFlann = o3d.KDTreeFlann
