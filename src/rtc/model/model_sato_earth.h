////////////////////////////////////////////////////////////////////////
// model_sato_earth.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_MODEL_SATO_EARTH_H
#define RTC_RAYTRACE_MODEL_SATO_EARTH_H

namespace rtc { namespace model { namespace plasma {

	/********************************************
	class sato_earth
	�@���̃��f���́A��������(2000�N����,PPARC)�ɂ����
	�ϑ����ꂽ�d�q���x���f�������ɍ쐬���Ă��܂��B

	*********************************************/
	class sato_earth : public basic_plasma_model
	{
	public:
		void test() const;

	protected:
		double getDensity( const vector& point ) const;

	public:
		struct SATO_PARAM
		{
			// ���̃p�����[�^�́A����(2000)�ɂ�����
			//�d�q���x�̓��o��
			//
			//   Ne = N0 * exp( -(z-z0)/(H0 + beta*z) );
			//
			// �ŗ��p�������̂ł���B
			// ���̍\���̂̃����o�ϐ����́A�㎮�̕ϐ����Ɉ�v����B
			// z �͒n�ォ��̋����ł���B
			//
			// ���̕������� getDensity()�ɃC���v�������g����Ă���B

			double
				MLT,  // ���C���[�J������[h]
				iLAT, // �s�ώ��C�ܓx  [deg]
				N0,   // �d�q���x�  [��/cm*3]
				z0,   // ���x�      [km]
				H0,
				beta;
		};

		class parameter : public SATO_PARAM
		{
		public:
			parameter(){ MLT = iLAT = N0 = z0 = H0 = beta = 0.0; };
			parameter( const SATO_PARAM& p );

			parameter& operator += ( const SATO_PARAM& p );
			parameter& operator =  ( const parameter& p );

			parameter operator +( const SATO_PARAM& p ) const;
			parameter operator *( double f ) const;
		};

	private:
		double getBaseDensity(
			const int MLT,
			const int iLAT,
			const double z
		) const;

		double getMLT ( const vector& point ) const;
		double getILAT( const vector& point ) const;

		parameter getBaseParam( int index ) const;

	private:
		// �p�����[�^�̊�ƂȂ�ϐ��Z�b�g.
		//9(MLT) x 40(iLAT)�̊�b�p�����[�^���i�[�����B
		//MLT,iLAT�Ɣz��C���f�b�N�X�̑Ή��́AgetTriBaseParamIndices()���Q�ƁB
		const static SATO_PARAM m_baseParam[];
		
		// �L���r�e�B�𐶐�����␳�I�u�W�F�N�g
		cavity m_cavity;
	};

}}}// namespace rtc

#endif//RTC_RAYTRACE_MODEL_SATO_EARTH_H
