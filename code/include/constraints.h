#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType { DISTANCE, COLLISION } ConstraintType;  // You can expand it for more constraints
typedef enum ConstraintEqualityType { EQUALITY, INEQUALITY } ConstraintEqualityType;

// there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for
// (x,y,z);
class Constraint {
   public:
    int m1, m2;  // Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
    int v1, v2;  // Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
    double invMass1, invMass2;      // inverse masses of two bodies
    double refValue;                // Reference values to use in the constraint, when needed (like distance)
    bool isUpper;                   // in case this is an inequality constraints, whether it's an upper or a lower bound
    RowVector3d refVector;          // Reference vector when needed (like vector)
    double CRCoeff;                 // velocity bias
    ConstraintType constraintType;  // The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change
                                    // after initialization!
    ConstraintEqualityType constraintEqualityType;  // whether the constraint is an equality or an inequality

    Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const bool _isUpper,
               const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2,
               const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff)
        : constraintType(_constraintType),
          constraintEqualityType(_constraintEqualityType),
          isUpper(_isUpper),
          m1(_m1),
          v1(_v1),
          m2(_m2),
          v2(_v2),
          invMass1(_invMass1),
          invMass2(_invMass2),
          refValue(_refValue),
          CRCoeff(_CRCoeff) {
        refVector = _refVector;
    }

    ~Constraint() {
    }

    // computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
    // The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
    // returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction
    // done)
    bool resolve_velocity_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions,
                                     const MatrixXd& currCOMVelocities, const MatrixXd& currAngVelocities,
                                     const Matrix3d& invInertiaTensor1, const Matrix3d& invInertiaTensor2,
                                     MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngVelocities, double tolerance) {
        /***************************TODO: implement this function**********************/

        correctedCOMVelocities = currCOMVelocities;
        correctedAngVelocities = currAngVelocities;

        const Vector3d p1 = currVertexPositions.row(0).transpose();
        const Vector3d p2 = currVertexPositions.row(1).transpose();
        const Vector3d J = (p1 - p2).normalized();
        const Vector3d arm1 = p1 - currCOMPositions.row(0).transpose();
        const Vector3d arm2 = p2 - currCOMPositions.row(1).transpose();
        const double norm =
            J.dot(((currCOMVelocities.row(0).transpose()) + ((Vector3d)currAngVelocities.row(0).transpose()).cross(arm1)) -
                  ((currCOMVelocities.row(1).transpose()) + ((Vector3d)currAngVelocities.row(1).transpose()).cross(arm2)));

        // Constraint is satisfied
        if (abs(norm) <= tolerance) return true;

        const double totalInvMass = invMass1 + invMass2;
        // Both objects are fixed
        if (totalInvMass == 0) return false;

        const Vector3d arm1XJ = arm1.cross(J);
        const Vector3d arm2XJ = arm2.cross(J);

        const double arm1Inertia = arm1XJ.transpose() * invInertiaTensor1 * arm1XJ;
        const double arm2Inertia = arm2XJ.transpose() * invInertiaTensor2 * arm2XJ;
        const double lambda = -(1 + CRCoeff) * norm / (totalInvMass + arm1Inertia + arm2Inertia);

        correctedCOMVelocities.row(0) += lambda * invMass1 * J.transpose();
        correctedCOMVelocities.row(1) -= lambda * invMass2 * J.transpose();

        correctedAngVelocities.row(0) += lambda * (invInertiaTensor1 * arm1XJ).transpose();
        correctedAngVelocities.row(1) -= lambda * (invInertiaTensor2 * arm2XJ).transpose();

        return false;
    }

    // projects the position unto the constraint
    // returns true if constraint was already good
    bool resolve_position_constraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions,
                                     MatrixXd& correctedCOMPositions, double tolerance) {
        /***************************TODO: implement this function**********************/

        correctedCOMPositions = currCOMPositions;

        const RowVector3d constraint = currConstPositions.row(1) - currConstPositions.row(0);
        const double dist = constraint.norm();
        const double diff = isUpper ? dist - refValue : refValue - dist;

        // Constraint is satisfied
        if (diff < tolerance) return true;

        const double totalInvMass = invMass1 + invMass2;
        // Both objects are fixed
        if (totalInvMass == 0) return false;

        const RowVector3d correction = (isUpper ? 1 : -1) * constraint / dist * diff / totalInvMass;

        correctedCOMPositions.row(0) += invMass1 * correction;
        correctedCOMPositions.row(1) -= invMass2 * correction;

        return false;
    }
};

#endif /* constraints_h */
