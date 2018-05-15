#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    float l = L/1.4142135623730951;
    float f1 = collThrustCmd;
    float f2 = momentCmd.x/l;
    float f3 = momentCmd.y/l;
    float f4 = momentCmd.z/kappa;
  
    cmd.desiredThrustsN[0] = 0.25 * ( f1 + f2 + f3 - f4 ); // front left
    cmd.desiredThrustsN[1] = 0.25 * ( f1 - f2 + f3 + f4 ); // front right
    cmd.desiredThrustsN[2] = 0.25 * ( f1 + f2 - f3 + f4 ); // rear left
    cmd.desiredThrustsN[3] = -0.25 * (-f1 + f2 + f3 + f4 ); // rear right

   //cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
   //cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
   //cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
   //cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    momentCmd.x = Ixx * kpPQR.x * ( pqrCmd.x - pqr.x ); 
    momentCmd.y = Iyy * kpPQR.y * ( pqrCmd.y - pqr.y ); 
    momentCmd.z = Izz * kpPQR.z * ( pqrCmd.z - pqr.z ); 
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float target_R13, target_R23, p_cmd, q_cmd;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  float c_d = collThrustCmd/mass;
  
  if ( collThrustCmd > 0.0 ) {
      target_R13 = -CONSTRAIN(accelCmd.x/c_d, -maxTiltAngle, maxTiltAngle); //#-min(max(accelCmd.x/c_d, -maxTiltAngle), maxTiltAngle)
      target_R23 = -CONSTRAIN(accelCmd.y/c_d, -maxTiltAngle, maxTiltAngle); //#-min(max(accelCmd.y/c_d, -maxTiltAngle), maxTiltAngle)
      
      p_cmd = (1/R(2, 2)) *
              (-R(1, 0) * kpBank * (R(0, 2)-target_R13) + 
               R(0, 0) * kpBank * (R(1, 2)-target_R23));
      q_cmd = (1/R(2, 2)) * 
              (-R(1, 1) * kpBank * (R(0, 2)-target_R13) + 
               R(0, 1) * kpBank * (R(1, 2)-target_R23));
  }
  else { //:  # Otherwise command no rate
      //print("negative thrust command")
      p_cmd = 0.0;
      q_cmd = 0.0;
  }
  
  pqrCmd.x = p_cmd;
  pqrCmd.y = q_cmd;
  pqrCmd.z = 0.0;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust, hdot_cmd, acceleration_cmd, tmp_thrust;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  hdot_cmd = kpPosZ * (posZCmd - posZ) + velZCmd;
  hdot_cmd = CONSTRAIN(hdot_cmd, -maxDescentRate, maxAscentRate);
  integratedAltitudeError += (posZCmd - posZ) * dt;
  acceleration_cmd = accelZCmd + kpVelZ*(hdot_cmd - velZ) + KiPosZ*integratedAltitudeError;
  thrust = -mass * ( acceleration_cmd - CONST_GRAVITY ) / R(2,2);
  //thrust = -(acceleration_cmd - CONST_GRAVITY) / R(2,2);
  tmp_thrust = thrust;

  //hdot_cmd = kpPosZ * (posZCmd - posZ);
  //hdot_cmd = CONSTRAIN(hdot_cmd, -maxAscentRate, maxDescentRate);
  //integratedAltitudeError += (posZCmd - posZ) * dt;
  //acceleration_cmd = accelZCmd + hdot_cmd + kpVelZ*(velZCmd - velZ) + KiPosZ*integratedAltitudeError;
  ////thrust = -(acceleration_cmd - CONST_GRAVITY) / R(2,2);
  //thrust = -mass * acceleration_cmd / R(2,2);

  //if (thrust > maxMotorThrust) { thrust = maxMotorThrust; }
  //else if (thrust < minMotorThrust) { thrust = minMotorThrust; }
  //printf ("tt=%f, t=%f, i=%f\n", tmp_thrust, thrust,  KiPosZ*integratedAltitudeError);
  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;
  V3F velocity_cmd;
  V3F accelCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  velocity_cmd = kpPosXY * ( posCmd - pos ) + velCmd;
  float velocity_norm = sqrt( velocity_cmd.x*velocity_cmd.x + velocity_cmd.y*velocity_cmd.y );
  if ( velocity_norm > maxSpeedXY ) {  velocity_cmd = velocity_cmd*maxSpeedXY/velocity_norm; }
  accelCmd = accelCmdFF + kpVelXY * ( velocity_cmd - vel );
  //accelCmd = kpVelXY * ( velCmd - vel ) + kpPosXY * ( posCmd - pos );
  accelCmd.z = 0;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  float pi=3.14159265;
  float yaw_err = yawCmd - yaw;
  yaw_err = fmodf( yaw_err, pi );
  yawRateCmd = kpYaw * yaw_err;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
