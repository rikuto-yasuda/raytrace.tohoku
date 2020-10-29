////////////////////////////////////////////////////////////////////////
// plasma_cavity.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_PLASMA_CAVITY_H
#define RTC_RAYTRACE_PLASMA_CAVITY_H

namespace rtc {
	
	// �v���Y�}�L���r�e�B���쐬���A�ݒu����B
	// ����ILAT���fp/fc���Œ�l���Ƃ�悤��
	// �v���Y�}���x������������B
	class cavity
	{
		friend class basic_plasma_model;

	public:
		cavity();
		cavity(
			double min_fpfc,    // cavity���S�ł́A��̂�fp/fc�l���w�肷��B
			double ilat_center, // cavity���S��ʂ�ILAT���w�肷��B
			double ilat_range,  // cavity�̈ܓx�����̔��l�����p�x�Ŏw�肷��B
			double mlt_center,  // cavity���S��ʂ�MLT���w�肷��B
			double mlt_range,   // cavity��LT�����̔��l���������Ŏw�肷��B
			double max_height,  // cavity���e�����y�ڂ��ō����x[km]���w�肷��B
			double bottom_height// cavity���S�̍��x[km]���w�肷��B
		);
		
		cavity& operator =( const cavity& c );
		double factor( const vector& ptr ) const;

		double getMinFpFc()                    const { return m_min_fpfc; }
		double getInvariantLatitudeCenterPtr() const { return m_ilat_center_ptr; }
		double getMagnetLocalTimeCenterPtr()   const { return m_mlt_center_ptr;  }
		double getInvariantLatitudeRange()     const { return m_ilat_range;      }
		double getMagnetLocalTimeRange()       const { return m_mlt_range;       }
		double getMaxHeight()                  const { return m_max_height;      }
		double getBottomHeight()               const { return m_bottomHeight;    }
		
		// �ȉ������I�Ɏg�p
	private:
		void create(
			const basic_planet& mother
		);
		
	private:
		double m_min_fpfc;
		double m_ilat_center_ptr;
		double m_mlt_center_ptr;
		double m_ilat_range;
		double m_mlt_range;
		double m_max_height;
		
		double m_bottomHeight;
		double m_scale;
		
		const basic_planet* m_mother;
	};
}// namespace rtc

#endif//RTC_RAYTRACE_PLASMA_CAVITY_H
