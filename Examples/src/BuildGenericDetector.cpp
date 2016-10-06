// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Examples/BuildGenericDetector.hpp"
#include <array>
#include <iostream>
#include <vector>
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Examples/GenericLayerBuilder.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Tools/LayerArrayCreator.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Tools/PassiveLayerBuilder.hpp"
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Tools/TrackingGeometryBuilder.hpp"
#include "ACTS/Tools/TrackingVolumeArrayCreator.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

std::unique_ptr<const Acts::TrackingGeometry>
buildGenericDetector(Logging::Level surfaceLLevel,
                     Logging::Level layerLLevel,
                     Logging::Level volumeLLevel,
                     size_t         stage)
{
  // configure surface array creator
  auto surfaceArrayCreator = std::make_shared<SurfaceArrayCreator>(
      getDefaultLogger("SurfaceArrayCreator", surfaceLLevel));
  // configure the layer creator that uses the surface array creator
  LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator            = std::make_shared<LayerCreator>(
      lcConfig, getDefaultLogger("LayerCreator", layerLLevel));
  // configure the layer array creator
  auto layerArrayCreator = std::make_shared<LayerArrayCreator>(
      getDefaultLogger("LayerArrayCreator", layerLLevel));
  // tracking volume array creator
  auto tVolumeArrayCreator = std::make_shared<TrackingVolumeArrayCreator>(
      getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator          = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper            = std::make_shared<CylinderVolumeHelper>(
      cvhConfig, getDefaultLogger("CylinderVolumeHelper", volumeLLevel));
  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<ITrackingVolumeBuilder>> volumeBuilders;

  //-------------------------------------------------------------------------------------
  // beam pipe
  //-------------------------------------------------------------------------------------
  // configure the beam pipe layer builder
  PassiveLayerBuilder::Config bplConfig;
  bplConfig.layerIdentification     = "BeamPipe";
  bplConfig.centralLayerRadii       = std::vector<double>(1, 19.);
  bplConfig.centralLayerHalflengthZ = std::vector<double>(1, 400.);
  bplConfig.centralLayerThickness   = std::vector<double>(1, 0.8);
  bplConfig.centralLayerMaterial = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
  auto beamPipeBuilder              = std::make_shared<PassiveLayerBuilder>(
      bplConfig, getDefaultLogger("BeamPipeLayerBuilder", layerLLevel));
  // create the volume for the beam pipe
  CylinderVolumeBuilder::Config bpvConfig;
  bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
  bpvConfig.volumeName           = "BeamPipe";
  bpvConfig.layerBuilder         = beamPipeBuilder;
  bpvConfig.layerEnvelopeR    = {1. * Acts::units::_mm, 1. * Acts::units::_mm};
  bpvConfig.buildToRadiusZero = true;
  bpvConfig.volumeSignature      = 0;
  auto beamPipeVolumeBuilder     = std::make_shared<CylinderVolumeBuilder>(
      bpvConfig, getDefaultLogger("BeamPipeVolumeBuilder", volumeLLevel));
  // add to the list of builders
  volumeBuilders.push_back(beamPipeVolumeBuilder);
  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // pixel detector
  //-------------------------------------------------------------------------------------
  // configure pixel layer builder
  GenericLayerBuilder::Config plbConfig;
  plbConfig.layerCreator        = layerCreator;
  plbConfig.layerIdentification = "Pixel";
  // fill necessary vectors for configuration
  //-------------------------------------------------------------------------------------
  // some prep work
  // envelope double
  std::pair<double, double> pcEnvelope(2., 2.);
  // Layer material properties - thickness, X0, L0, A, Z, Rho
  MaterialProperties pcmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
  // Module material - X0, L0, A, Z, Rho
  Material pcMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

  // standard, an approach envelope
  plbConfig.approachSurfaceEnvelope = 0.5;

  // STAGE 0 --- 1 pixel layer detector for debugging
  if (stage == 0) {
    // configure the central barrel
    plbConfig.centralLayerBinMultipliers = {1, 1};
    plbConfig.centralLayerRadii                 = {29.};
    plbConfig.centralLayerEnvelopes             = {pcEnvelope};
    plbConfig.centralLayerMaterialConcentration = {1};
    plbConfig.centralLayerMaterialProperties    = {pcmProperties};
    plbConfig.centralModuleBinningSchema        = {{16, 13}};
    plbConfig.centralModuleTiltPhi              = {0.18};
    plbConfig.centralModuleHalfX                = {8.4};
    plbConfig.centralModuleHalfY                = {32.};
    plbConfig.centralModuleThickness            = {0.15};
    plbConfig.centralModuleMaterial             = {pcMaterial};
  } else {
    // configure the central barrel
    plbConfig.centralLayerBinMultipliers = {1, 1};
    // STAGE > 0 --- 4 pixel layers
    plbConfig.centralLayerRadii = {29., 55., 88., 124.};
    plbConfig.centralLayerEnvelopes
        = {pcEnvelope, pcEnvelope, pcEnvelope, pcEnvelope};
    plbConfig.centralLayerMaterialConcentration = {1, 1, 1, 1};
    plbConfig.centralLayerMaterialProperties
        = {pcmProperties, pcmProperties, pcmProperties, pcmProperties};
    plbConfig.centralModuleBinningSchema
        = {{16, 13}, {24, 13}, {38, 13}, {60, 13}};
    plbConfig.centralModuleTiltPhi   = {0.18, 0.18, 0.2, 0.2};
    plbConfig.centralModuleHalfX     = {8.4, 8.4, 8.4, 8.4};
    plbConfig.centralModuleHalfY     = {32., 32., 32., 32.};
    plbConfig.centralModuleThickness = {0.15, 0.15, 0.15, 0.15};
    plbConfig.centralModuleMaterial
        = {pcMaterial, pcMaterial, pcMaterial, pcMaterial};
  }
  // no frontside/backside
  plbConfig.centralModuleFrontsideStereo = {};
  plbConfig.centralModuleBacksideStereo  = {};
  plbConfig.centralModuleBacksideGap     = {};
  // mPositions
  std::vector<std::vector<Vector3D>> centralModulePositions;
  for (size_t plb = 0; plb < plbConfig.centralLayerRadii.size(); ++plb) {
    // call the helper function
    centralModulePositions.push_back(
        modulePositionsCylinder(plbConfig.centralLayerRadii[plb],
                                0.5,  // 1 mm stagger
                                plbConfig.centralModuleHalfY[plb],
                                2.,  // 2 mm module overlap
                                plbConfig.centralModuleBinningSchema[plb]));
  }
  plbConfig.centralModulePositions = centralModulePositions;
  // STAGE = 0 - 1 disk detector
  if (stage == 0) {
    // configure the endcaps
    plbConfig.posnegLayerBinMultipliers        = {1, 1};
    plbConfig.posnegLayerPositionsZ            = {500};
    plbConfig.posnegLayerEnvelopeR             = {1.};
    plbConfig.posnegLayerMaterialConcentration = {1};
    plbConfig.posnegLayerMaterialProperties    = {pcmProperties};
    plbConfig.posnegModuleMinHalfX             = {{8.4}};
    plbConfig.posnegModuleMaxHalfX             = {};
    plbConfig.posnegModuleHalfY                = {{42.}};
    plbConfig.posnegModulePhiBins              = {{48}};
    plbConfig.posnegModuleThickness            = {{0.15}};
    plbConfig.posnegModuleMaterial             = {{pcMaterial}};
  } else {
    // STAGE > 0 - pixel endcap detector
    // configure the endcaps
    plbConfig.posnegLayerBinMultipliers        = {1, 1};
    plbConfig.posnegLayerPositionsZ            = {500., 580., 680., 800.};
    plbConfig.posnegLayerEnvelopeR             = {1., 1., 1., 1.};
    plbConfig.posnegLayerMaterialConcentration = {1, 1, 1, 1};
    plbConfig.posnegLayerMaterialProperties
        = {pcmProperties, pcmProperties, pcmProperties, pcmProperties};
    plbConfig.posnegModuleMinHalfX  = {{8.4}, {8.4}, {8.4}, {8.4}};
    plbConfig.posnegModuleMaxHalfX  = {};
    plbConfig.posnegModuleHalfY     = {{42.}, {42.}, {42.}, {42.}};
    plbConfig.posnegModulePhiBins   = {{48}, {48}, {48}, {48}};
    plbConfig.posnegModuleThickness = {{0.15}, {0.15}, {0.15}, {0.15}};
    plbConfig.posnegModuleMaterial
        = {{pcMaterial}, {pcMaterial}, {pcMaterial}, {pcMaterial}};
  }
  // no frontside/backside
  plbConfig.posnegModuleFrontsideStereo = {};
  plbConfig.posnegModuleBacksideStereo  = {};
  plbConfig.posnegModuleBacksideGap     = {};
  // mPositions
  std::vector<std::vector<std::vector<Vector3D>>> posnegModulePositions;
  for (size_t id = 0; id < plbConfig.posnegLayerPositionsZ.size(); ++id) {
    posnegModulePositions.push_back(
        modulePositionsDisc(plbConfig.posnegLayerPositionsZ[id],
                            2.0,
                            0.5,
                            29.,
                            120.,
                            plbConfig.posnegModulePhiBins[id],
                            plbConfig.posnegModuleHalfY[id]));
  }
  plbConfig.posnegModulePositions = posnegModulePositions;

  // define the builder
  auto pixelLayerBuilder = std::make_shared<GenericLayerBuilder>(
      plbConfig, getDefaultLogger("PixelLayerBuilder", layerLLevel));
  //-------------------------------------------------------------------------------------
  // build the pixel volume
  CylinderVolumeBuilder::Config pvbConfig;
  pvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
  pvbConfig.volumeName           = "Pixel";
  pvbConfig.buildToRadiusZero    = false;
  pvbConfig.layerEnvelopeR = {1. * Acts::units::_mm, 5. * Acts::units::_mm};
  pvbConfig.layerBuilder         = pixelLayerBuilder;
  pvbConfig.volumeSignature      = 0;
  auto pixelVolumeBuilder        = std::make_shared<CylinderVolumeBuilder>(
      pvbConfig, getDefaultLogger("PixelVolumeBuilder", volumeLLevel));
  // add to the list of builders
  volumeBuilders.push_back(pixelVolumeBuilder);

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // SHORT strip detector
  //-------------------------------------------------------------------------------------
  // STAGE > 1 - short strip detector added
  if (stage > 1) {
    // first add a Pixel Support Tube
    //-------------------------------------------------------------------------------------
    // Pixel Support Tybe (PST)
    //-------------------------------------------------------------------------------------
    PassiveLayerBuilder::Config pstConfig;
    pstConfig.layerIdentification     = "PST";
    pstConfig.centralLayerRadii       = std::vector<double>(1, 175.);
    pstConfig.centralLayerHalflengthZ = std::vector<double>(1, 1200.);
    pstConfig.centralLayerThickness   = std::vector<double>(1, 1.8);
    pstConfig.centralLayerMaterial
        = {Material(352.8, 407., 9.012, 4., 1.848e-3)};
    auto pstBuilder = std::make_shared<PassiveLayerBuilder>(
        pstConfig, getDefaultLogger("PstBuilder", layerLLevel));
    // create the volume for the beam pipe
    CylinderVolumeBuilder::Config pstvolConfig;
    pstvolConfig.trackingVolumeHelper = cylinderVolumeHelper;
    pstvolConfig.volumeName           = "PST";
    pstvolConfig.buildToRadiusZero    = false;
    pstvolConfig.layerBuilder         = pstBuilder;
    pstvolConfig.volumeSignature      = 0;
    auto pstVolumeBuilder             = std::make_shared<CylinderVolumeBuilder>(
        pstvolConfig, getDefaultLogger("PstVolumeBuilder", volumeLLevel));
    // add to the detector builds
    // @TODO check why this is not yet working
    // volumeBuilders.push_back(pstVolumeBuilder);

    // STRIPS
    // ----------------------------------------------------------------------------
    // configure short strip layer builder
    GenericLayerBuilder::Config sslbConfig;
    sslbConfig.layerCreator        = layerCreator;
    sslbConfig.layerIdentification = "SStrip";
    // fill necessary vectors for configuration
    //-------------------------------------------------------------------------------------
    // some prep work
    // envelope double
    std::pair<double, double> ssEnvelope(2., 2.);
    // Layer material properties - thickness, X0, L0, A, Z, Rho
    MaterialProperties ssmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
    // Module material - X0, L0, A, Z, Rho
    Material ssMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

    // configure the central barrel
    sslbConfig.centralLayerBinMultipliers = {1, 1};
    sslbConfig.centralLayerRadii          = {220., 350., 500.};
    sslbConfig.centralLayerEnvelopes = {ssEnvelope, ssEnvelope, ssEnvelope};
    sslbConfig.centralLayerMaterialConcentration = {1, 1, 1};
    sslbConfig.centralLayerMaterialProperties
        = {ssmProperties, ssmProperties, ssmProperties};
    sslbConfig.centralModuleBinningSchema = {{42, 12}, {64, 12}, {84, 12}};
    sslbConfig.centralModuleTiltPhi       = {-0.15, -0.15, -0.15};
    sslbConfig.centralModuleHalfX         = {18.2, 18.2, 18.2};
    sslbConfig.centralModuleHalfY         = {68., 68., 68.};
    sslbConfig.centralModuleThickness     = {0.25, 0.25, 0.25};
    sslbConfig.centralModuleMaterial
        = {ssMaterial, ssMaterial, ssMaterial, ssMaterial};
    sslbConfig.centralModuleFrontsideStereo = {-0.02, -0.02, -0.02};
    sslbConfig.centralModuleBacksideStereo  = { 0.02, 0.02, 0.02};
    sslbConfig.centralModuleBacksideGap     = {2., 2., 2.};
    // mPositions
    std::vector<std::vector<Vector3D>> centralModulePositions;
    for (size_t sslb = 0; sslb < sslbConfig.centralLayerRadii.size(); ++sslb) {
      // call the helper function
      centralModulePositions.push_back(
          modulePositionsCylinder(sslbConfig.centralLayerRadii[sslb],
                                  0.5,  // 1 mm stagger
                                  sslbConfig.centralModuleHalfY[sslb],
                                  2.,  // 2 mm module overlap
                                  sslbConfig.centralModuleBinningSchema[sslb]));
    }
    sslbConfig.centralModulePositions = centralModulePositions;

    // configure the endcaps
    std::vector<double>   mrMinHx    = {16.4, 24.2, 32.2};
    std::vector<double>   mrMaxHx    = {24.2, 32.2, 40.0};
    std::vector<double>   mrHy       = {48., 48., 48.};
    std::vector<int>      mPhiBins   = {42, 58, 72};
    std::vector<double>   mThickness = {0.2, 0.2, 0.2};
    std::vector<Material> mMaterial  = {ssMaterial, ssMaterial, ssMaterial};
    std::vector<double>   mfStereo   = {-0.02, -0.02, -0.02};
    std::vector<double>   mbStereo   = {0.02, 0.02, 0.02};
    std::vector<double>   mfbGap     = {2., 2., 2.};

    sslbConfig.posnegLayerBinMultipliers        = { 1, 2 };
    sslbConfig.posnegLayerPositionsZ
        = {880., 1100., 1300., 1550., 1800., 2200.};
    size_t nposnegs                 = sslbConfig.posnegLayerPositionsZ.size();
    sslbConfig.posnegLayerEnvelopeR = std::vector<double>(nposnegs, 5.);
    sslbConfig.posnegLayerMaterialConcentration = std::vector<int>(nposnegs, 1);
    sslbConfig.posnegLayerMaterialProperties
        = std::vector<MaterialProperties>(nposnegs, ssmProperties);
    sslbConfig.posnegModuleMinHalfX
        = std::vector<std::vector<double>>(nposnegs, mrMinHx);
    sslbConfig.posnegModuleMaxHalfX
        = std::vector<std::vector<double>>(nposnegs, mrMaxHx);
    sslbConfig.posnegModuleHalfY
        = std::vector<std::vector<double>>(nposnegs, mrHy);
    sslbConfig.posnegModulePhiBins
        = std::vector<std::vector<int>>(nposnegs, mPhiBins);
    sslbConfig.posnegModuleThickness
        = std::vector<std::vector<double>>(nposnegs, mThickness);
    sslbConfig.posnegModuleMaterial
        = std::vector<std::vector<Material>>(nposnegs, mMaterial);
    sslbConfig.posnegModuleFrontsideStereo
        = std::vector<std::vector<double>>(nposnegs, mfStereo);
    sslbConfig.posnegModuleBacksideStereo
        = std::vector<std::vector<double>>(nposnegs, mbStereo);
    sslbConfig.posnegModuleBacksideGap
        = std::vector<std::vector<double>>(nposnegs, mfbGap);
    // mPositions
    std::vector<std::vector<std::vector<Vector3D>>> posnegModulePositions;
    for (size_t id = 0; id < sslbConfig.posnegLayerPositionsZ.size(); ++id) {
      posnegModulePositions.push_back(
          modulePositionsDisc(sslbConfig.posnegLayerPositionsZ[id],
                              2.0,
                              0.5,
                              220.,
                              500.,
                              sslbConfig.posnegModulePhiBins[id],
                              sslbConfig.posnegModuleHalfY[id]));
    }
    sslbConfig.posnegModulePositions = posnegModulePositions;

    // define the builder
    auto sstripLayerBuilder = std::make_shared<GenericLayerBuilder>(
        sslbConfig, getDefaultLogger("SStripLayerBuilder", layerLLevel));
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    CylinderVolumeBuilder::Config ssvbConfig;
    ssvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    ssvbConfig.volumeName           = "SStrip";
    ssvbConfig.buildToRadiusZero    = false;
    ssvbConfig.layerBuilder         = sstripLayerBuilder;
    ssvbConfig.volumeSignature      = 0;
    auto sstripVolumeBuilder        = std::make_shared<CylinderVolumeBuilder>(
        ssvbConfig, getDefaultLogger("SStripVolumeBuilder", volumeLLevel));

    //-------------------------------------------------------------------------------------
    // add to the list of builders
    volumeBuilders.push_back(sstripVolumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // LONG strip detector
  //-------------------------------------------------------------------------------------
  // STAGE > 2 - long strip detector added
  if (stage > 2) {
    // configure short strip layer builder
    GenericLayerBuilder::Config lslbConfig;
    lslbConfig.layerCreator        = layerCreator;
    lslbConfig.layerIdentification = "LStrip";
    // fill necessary vectors for configuration
    //-------------------------------------------------------------------------------------
    // some prep work
    // envelope double
    std::pair<double, double> lsEnvelope(2., 2.);
    // Layer material properties - thickness, X0, L0, A, Z, Rho
    MaterialProperties lsmProperties(1., 95.7, 465.2, 28.03, 14., 2.32e-3);
    // Module material - X0, L0, A, Z, Rho
    Material lsMaterial(95.7, 465.2, 28.03, 14., 2.32e-3);

    // configure the central barrel
    lslbConfig.centralLayerBinMultipliers        = {1, 1};
    lslbConfig.centralLayerRadii                 = {680., 980.};
    lslbConfig.centralLayerEnvelopes             = {lsEnvelope, lsEnvelope};
    lslbConfig.centralLayerMaterialConcentration = {1, 1};
    lslbConfig.centralLayerMaterialProperties = {lsmProperties, lsmProperties};
    lslbConfig.centralModuleBinningSchema     = {{64, 16}, {98, 16}};
    lslbConfig.centralModuleTiltPhi           = {-0.15, -0.15};
    lslbConfig.centralModuleHalfX             = {42., 42.};
    lslbConfig.centralModuleHalfY             = {76., 76.};
    lslbConfig.centralModuleThickness         = {0.25, 0.25, 0.25};
    lslbConfig.centralModuleMaterial
        = {lsMaterial, lsMaterial, lsMaterial, lsMaterial};
    lslbConfig.centralModuleFrontsideStereo = {-0.02, -0.02};
    lslbConfig.centralModuleBacksideStereo  = {0.02, 0.02};
    lslbConfig.centralModuleBacksideGap     = {2., 2.};
    // mPositions
    std::vector<std::vector<Vector3D>> centralModulePositions;
    for (size_t lslb = 0; lslb < lslbConfig.centralLayerRadii.size(); ++lslb) {
      // call the helper function
      centralModulePositions.push_back(
          modulePositionsCylinder(lslbConfig.centralLayerRadii[lslb],
                                  0.5,  // 1 mm stagger
                                  lslbConfig.centralModuleHalfY[lslb],
                                  2.,  // 2 mm module overlap
                                  lslbConfig.centralModuleBinningSchema[lslb]));
    }
    lslbConfig.centralModulePositions = centralModulePositions;

    // configure the endcaps
    std::vector<double>   mrMinHx    = { 42., 42., 42. } ;
    std::vector<double>   mrMaxHx    = { 56., 56., 56.} ;
    std::vector<double>   mrHy       = { 64., 64., 64. } ;
    std::vector<int>      mPhiBins   = { 64, 78, 98 };
    std::vector<double>   mThickness = { 0.25, 0.25, 0.25  };
    std::vector<Material> mMaterial  = { lsMaterial, lsMaterial, lsMaterial };
    std::vector<double>   mfStereo   = { -0.02, -0.02, -0.02 };
    std::vector<double>   mbStereo   = { 0.02, 0.02, 0.02 };
    std::vector<double>   mfbGap     = { 2., 2., 2. };

    // endcap
    lslbConfig.posnegLayerBinMultipliers        = { 1, 2 };
    lslbConfig.posnegLayerPositionsZ            = { 1380., 1680., 2180. };
    size_t nposnegs =  lslbConfig.posnegLayerPositionsZ.size();
    lslbConfig.posnegLayerEnvelopeR             =
    std::vector<double>(nposnegs,5.);
    lslbConfig.posnegLayerMaterialConcentration = std::vector<int>(nposnegs,1);
    lslbConfig.posnegLayerMaterialProperties    =
    std::vector<MaterialProperties>(nposnegs,lsmProperties);
    lslbConfig.posnegModuleMinHalfX             = std::vector<
    std::vector<double> >(nposnegs, mrMinHx);
    lslbConfig.posnegModuleMaxHalfX             = std::vector<
    std::vector<double> >(nposnegs, mrMaxHx);
    lslbConfig.posnegModuleHalfY                = std::vector<
    std::vector<double> >(nposnegs, mrHy);
    lslbConfig.posnegModulePhiBins              = std::vector< std::vector< int
    > > (nposnegs, mPhiBins);
    lslbConfig.posnegModuleThickness            = std::vector<
    std::vector<double> >(nposnegs, mThickness);
    lslbConfig.posnegModuleMaterial             = std::vector<
    std::vector<Material> >(nposnegs, mMaterial);
    lslbConfig.posnegModuleFrontsideStereo      = std::vector<
    std::vector<double> >(nposnegs, mfStereo);
    lslbConfig.posnegModuleBacksideStereo       = std::vector<
    std::vector<double> >(nposnegs, mbStereo);
    lslbConfig.posnegModuleBacksideGap          = std::vector<
    std::vector<double> >(nposnegs, mfbGap);
    // mPositions
    std::vector<std::vector<std::vector<Vector3D>>> posnegModulePositions;
    for (size_t id = 0; id < lslbConfig.posnegLayerPositionsZ.size(); ++id) {
      posnegModulePositions.push_back(
          modulePositionsDisc(lslbConfig.posnegLayerPositionsZ[id],
                              2.0,
                              0.5,
                              600.,
                              980.,
                              lslbConfig.posnegModulePhiBins[id],
                              lslbConfig.posnegModuleHalfY[id]));
    }
    lslbConfig.posnegModulePositions = posnegModulePositions;

    // define the builder
    auto lstripLayerBuilder = std::make_shared<GenericLayerBuilder>(
        lslbConfig, getDefaultLogger("LStripLayerBuilder", layerLLevel));
    //-------------------------------------------------------------------------------------
    // build the pixel volume
    CylinderVolumeBuilder::Config lsvbConfig;
    lsvbConfig.trackingVolumeHelper = cylinderVolumeHelper;
    lsvbConfig.volumeName           = "LStrip";
    lsvbConfig.buildToRadiusZero    = false;
    lsvbConfig.layerBuilder         = lstripLayerBuilder;
    lsvbConfig.volumeSignature      = 0;
    auto lstripVolumeBuilder        = std::make_shared<CylinderVolumeBuilder>(
        lsvbConfig, getDefaultLogger("LStripVolumeBuilder", volumeLLevel));
    // add to the list of builders
    volumeBuilders.push_back(lstripVolumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  TrackingGeometryBuilder::Config tgConfig;
  tgConfig.trackingVolumeBuilders = volumeBuilders;
  tgConfig.trackingVolumeHelper   = cylinderVolumeHelper;
  auto cylinderGeometryBuilder
      = std::make_shared<const TrackingGeometryBuilder>(
          tgConfig, getDefaultLogger("TrackerGeometryBuilder", volumeLLevel));
  return cylinderGeometryBuilder->trackingGeometry();
}

/// helper method for cylinder
std::vector<Acts::Vector3D>
modulePositionsCylinder(double radius,
                        double zStagger,
                        double moduleHalfLength,
                        double lOverlap,
                        const std::pair<int, int>& binningSchema)
{
  int nPhiBins = binningSchema.first;
  int nZbins   = binningSchema.second;
  // prepare the return value
  std::vector<Vector3D> mPositions;
  mPositions.reserve(nPhiBins * nZbins);
  // prep work
  double phiStep = 2 * M_PI / (nPhiBins);
  double minPhi  = -M_PI + 0.5 * phiStep;
  double zStart  = -0.5 * (nZbins - 1) * (2 * moduleHalfLength - lOverlap);
  double zStep   = 2 * fabs(zStart) / (nZbins - 1);
  // loop over the bins
  for (size_t zBin = 0; zBin < size_t(nZbins); ++zBin) {
    // prepare z and r
    double moduleZ = zStart + zBin * zStep;
    double moduleR
        = (zBin % 2) ? radius - 0.5 * zStagger : radius + 0.5 * zStagger;
    for (size_t phiBin = 0; phiBin < size_t(nPhiBins); ++phiBin) {
      // calculate the current phi value
      double modulePhi = minPhi + phiBin * phiStep;
      mPositions.push_back(Vector3D(
          moduleR * cos(modulePhi), moduleR * sin(modulePhi), moduleZ));
    }
  }
  return mPositions;
}

/// helper method for disc
std::vector<std::vector<Acts::Vector3D>>
modulePositionsDisc(double                     z,
                    double                     ringStagger,
                    double                     phiStagger,
                    double                     innerRadius,
                    double                     outerRadius,
                    const std::vector<int>&    discBinning,
                    const std::vector<double>& moduleHalfLength)
{
  // calculate the radii
  std::vector<double> radii;
  // calculate the radial borders
  std::vector<double> radialBoarders;
  // the radial span of the disc
  double deltaR = outerRadius - innerRadius;
  // quick exits
  if (discBinning.size() == 1) {
    radii.push_back(0.5 * (innerRadius + outerRadius));
    radialBoarders = {innerRadius, outerRadius};
  } else {
    double totalLength = 0;
    // sum up the total length
    for (auto& mhlength : moduleHalfLength) totalLength += 2 * mhlength;
    // now calculate the overlap (equal pay)
    double rOverlap = (totalLength - deltaR) / (moduleHalfLength.size() - 1);
    // and now fill the radii and gaps
    double lastR  = innerRadius;
    double lastHl = 0.;
    double lastOl = 0.;
    // remember the radial boarders
    radialBoarders.push_back(innerRadius);
    // now calculate
    for (auto& mhlength : moduleHalfLength) {
      // calculate the radius
      radii.push_back(lastR + lastHl - lastOl + mhlength);
      lastR  = radii[radii.size() - 1];
      lastOl = rOverlap;
      lastHl = mhlength;
      // and register the radial boarder
      radialBoarders.push_back(lastR + 2 * lastHl - 0.5 * lastOl);
    }
  }
  // now prepare the return method
  std::vector<std::vector<Vector3D>> mPositions;
  for (size_t ir = 0; ir < radii.size(); ++ir) {
    // generate the z value
    double rz = radii.size() == 1 ? z : (ir % 2 ? z - 0.5 * ringStagger
                                                : z + 0.5 * ringStagger);
    // fill the ring positions
    mPositions.push_back(
        modulePositionsRing(rz, radii[ir], phiStagger, discBinning[ir]));
  }
  return mPositions;
}

/// Helper method for positioning
std::vector<Acts::Vector3D>
modulePositionsRing(double z, double radius, double phiStagger, int nPhiBins)
{
  // create and fill the positions
  std::vector<Vector3D> rPositions;
  rPositions.reserve(nPhiBins);
  // prep work
  double phiStep = 2 * M_PI / (nPhiBins);
  double minPhi  = -M_PI + 0.5 * phiStep;
  // phi loop
  for (size_t iphi = 0; iphi < size_t(nPhiBins); ++iphi) {
    double phi = minPhi + iphi * phiStep;
    double rz  = iphi % 2 ? z - 0.5 * phiStagger : z + 0.5 * phiStagger;
    rPositions.push_back(Vector3D(radius * cos(phi), radius * sin(phi), rz));
  }
  return rPositions;
}

}  // end of namespace Acts