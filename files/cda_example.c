// vertex parameters
double vertex[3];
double cda;
double midpoint[3];          // track midpoint coordinates at z = (vertex[2] + DCHbz)/2

// beam parameters
double beamCorrX, beamCorrY;          // beam slopes
double beamOffsetX, beamOffsetY;      // beam offsets at z=0

int          ntrack = sevt->Ntrack;

// Now proceed with the normal "One Good Track" identification like before v55
for (i=0; i<ntrack; i++)
  {    
    if (ntrack == 1)
      break;
    
    /////// v63: Get Beam parameters depending on beam type 
    if (sevt->track[i].q < 0)       // negative track == K- beam
      {	
	beamCorrX = abcog_params.pkdxdzm;
	beamOffsetX = abcog_params.pkxoffm;
	beamCorrY = abcog_params.pkdydzm;
	beamOffsetY = abcog_params.pkyoffm;
      } 	  
    else                            // positive track == K+ beam 
      {
	beamCorrX = abcog_params.pkdxdzp;
	beamOffsetX = abcog_params.pkxoffp;
	beamCorrY = abcog_params.pkdydzp;
	beamOffsetY = abcog_params.pkyoffp;
      }
    
    
    ///////////////////////////////////////////////////////////////////////////////////
    /////// Vertex-Berechnung ohne Bluefield-Korrektur
    ///////////////////////////////////////////////////////////////////////////////////
    
    // First define the coordinates of track and beam at DCH1 + the vectors, needed as input for vertex calculation
    
    // Track coordinates (p1 + v1)
    bdxdz_track = sevt->track[i].bdxdz;
    bdydz_track = sevt->track[i].bdydz;
    bdzdz_track = 1.;
    bx_track = sevt->track[i].bx;         // x coordinate before magnet (genauer: direkt (ca. 8cm) vor DCH1)
    by_track = sevt->track[i].by;         // y coordinate before magnet (genauer: direkt (ca. 8cm) vor DCH1)
    bz_track = DCHbz;
    
    v1[0] = bdxdz_track;
    v1[1] = bdydz_track;
    v1[2] = bdzdz_track;
    p1[0] = bx_track;
    p1[1] = by_track;
    p1[2] = bz_track;
    
    // Beam coordinates (p2 + v2)
    // Kaon mit Flugrichtung entlang Strahlachse 
    bdxdz_kaon = beamCorrX;
    bdydz_kaon = beamCorrY;
    bdzdz_kaon = 1.;
    bx_kaon = beamOffsetX + beamCorrX*DCHbz;
    by_kaon = beamOffsetY + beamCorrY*DCHbz;
    bz_kaon = DCHbz;
    
    v2[0] = bdxdz_kaon;
    v2[1] = bdydz_kaon;
    v2[2] = bdzdz_kaon;
    p2[0] = bx_kaon;
    p2[1] = by_kaon;
    p2[2] = bz_kaon;
    
    // Call routine for vertex calculation
    closap_double_(p1, p2, v1, v2, &cda, vertex);
  }
