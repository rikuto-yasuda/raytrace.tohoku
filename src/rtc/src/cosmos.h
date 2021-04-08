////////////////////////////////////////////////////////////////////////
// cosmos.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_COSMOS_H
#define RTC_RAYTRACE_COSMOS_H

namespace rtc {

	// cosmos ------------------------------------
	// �F���𒊏ۉ�����B
	// �F���̈ӎu�ɂ���Č������ݏo�����B
	// �F���ɂ͘f�������݂���B�����̘f����
	// �������ꂽ��A�F���Ɂu�o�^�v����K�v������B
	// �����v���Y�}���x�Ȃǂ́A�F���ɖ₢���킹�邱�ƂŎ擾�ł���B
	// ����ɂ���āA�s�v�ł͂���Ǝv�����A�����I��
	// �����̘f���A���邢�͉q���ɂ�镡���I�ȃ��f���̉��Z���\�ƂȂ�B
	// ���A���݂͓o�^�ł���f���͈�݂̂ł���B
	//
	// libraytrace�ł̓V���O���g���Ƃ��đ��݂��Ȃ���΂Ȃ�Ȃ��B
	//
	class cosmos
	{
	public:
		typedef std::set<ray*>    rays_t;
		
	public:
		cosmos(
			const int year,
			const int month,
			const int mday,
			const int hour,
			const int min,
			const int sec = 0 // �F���̎����� UT �Ŏw�肷��B
		);
		virtual ~cosmos();
		
		// ���� ------------------------------------------------------
		ray* createRay( const wave_parameter& wparam );
		ray* createRay(
			wave_parameter::wave_mode mode,   // �g���̃��[�h�BLO_MODE��RX_MODE���w�肷��B
			double                    freq,   // �g���̎��g��[Hz]���w�肷��B
			double              prec  = 3.74e-4, // step�O��̔�̋��e�����w�肷��B
			double              lstep = cnst::c, // ������1step�Ői�ލő咷���w��
			const double timeStep_max = 1e0,  // 1step�Ői�ގ��Ԃ̍ő�l���w��
			const double timeStep_min = 1e-54 // 1step�Ői�ގ��Ԃ̍ŏ��l���w��
		);
		
		// �������ꂽ�����ւ̃A�N�Z�X
		const rays_t& getRays() const
		{ return m_rays; }
		
		// �����̔j��
		void eraseRay( ray* pray );


		// ���� ------------------------------------------------------
		const std::tm& getUniversalTime() const
		{ return m_UniversalTime; }

		void setUniversalTime( const std::tm& ut )
		{ std::memcpy( &m_UniversalTime, &ut, sizeof(ut) ); }


		// �f�� ------------------------------------------------------
		const basic_planet& getPlanet() const
		{ return *m_planet; }
		
		bool registerPlanet( basic_planet& planet );


		// �����I���l�ւ̃A�N�Z�X ------------------------------------
		// �ʒur�ł̎���x�N�g����Ԃ��B
		vector getMagnetField( const vector& r ) const;
		
		// ������z���s��Ƃ��ĕԂ��B
		// �l�͈ȉ��̍s��ŕԂ��悤�ɂ��Ȃ���΂Ȃ�Ȃ��B
		// | dBx/dx dBx/dy dBx/dz |
		// | dBy/dx dBy/dy dBy/dz |
		// | dBz/dx dBz/dy dBz/dz |
		//
		// cosmos�ł́Acosmos::getMagnetField()�𕡐���Ăяo���A
		// ����̌X���𐔒l�I�Ɍv�Z���郁�\�b�h���L�q����Ă���B
		// ��͓I�ɋ��߂����̂ł͂Ȃ��̂Œ��ӂ��邱�ƁB
		matrix getDerivativeB( const vector& r ) const;
		
		// �v���Y�}���x�̌��z���x�N�g���ŕԂ��B
		// cosmos�ł́Acosmos::getPlasmaDensity()�𕡐���Ăяo���A
		// ���x���z�𐔒l�I�Ɍv�Z���郁�\�b�h���L�q����Ă���B
		// ��͓I�ɋ��߂����̂ł͂Ȃ��̂Œ��ӂ��邱�ƁB
		// ���z�̊Ԋu�́ARTC_DERIVATIVE_DISTANCE �Œ�`�����B
		vector getDerivativeDensity( const vector& r ) const;
		
		// �ʒur�ł̃v���Y�}���x��Ԃ��B
		double getPlasmaDensity( const vector& r ) const;
		
		// �ʒur�ł̕W���B
		double getHight( const vector& r ) const;
		
		// �ʒur�ł̃v���Y�}�p���g���̓��l��Hight�Ԃ��B
		double getSquaredPlasmaFreq( const vector& r ) const;

		// �ʒur�ł̃T�C�N���g�����p���g���̓��l��Ԃ��B
		double getSquaredCycloFreq ( const vector& r ) const;

	private:
		std::tm             m_UniversalTime;

		rays_t              m_rays;
		basic_planet*       m_planet;

		// �}���`�X���b�h�Ή�
		boost::mutex m_lock;
	};

}
#endif // !defined(AFX_ENV_H__933C5528_1711_4028_A437_6E3D8421B3E3__INCLUDED_)
