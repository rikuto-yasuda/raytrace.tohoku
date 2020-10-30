////////////////////////////////////////////////////////////////////////
// earth.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_PLANET_EARTH_H
#define RTC_RAYTRACE_PLANET_EARTH_H

namespace rtc { namespace planet {
	
	// �n����\���f���N���X
	class earth : public basic_planet
	{
	public:
		earth(
			basic_magnet_model&   mag,        // ���ꃂ�f���̃C���X�^���X���w�肷��B
			basic_plasma_model&  plsm         // �v���Y�}���f���̃C���X�^���X���w�肷��B
		);

		// ��]�s�� --------------------------------------------------
		virtual matrix getGEI2GEO() const;
		virtual matrix getGEI2GSE() const;
		virtual matrix getGSE2GSM() const;
		virtual matrix getGSM2SM () const;

	protected:
		double getMJD() const;
	};
	
}}// namespace

#endif//RTC_RAYTRACE_PLANET_EARTH_H
