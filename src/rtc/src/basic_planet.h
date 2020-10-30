////////////////////////////////////////////////////////////////////////
// basic_planet.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_BASIC_PLANET_H
#define RTC_RAYTRACE_BASIC_PLANET_H

namespace rtc {
	
	// basic_planet ------------------------------
	// �f���̊�{�N���X�B
	// �f���͎���ƃv���Y�}�̊e���f���ɑ΂�
	// has a �̊֌W�ɂ���A�ǂ̃��f�����g�p���邩��
	// �R���X�g���N�^�Ŏw�肵�Ȃ���΂Ȃ�Ȃ��B
	//
	class basic_planet
	{
		friend class cosmos;

	public:
		// �n���A����ю����̏����i�[����B
		class axis_info
		{
		public:
			axis_info(
				const double magnet_latitude, // �����̈ܓx [deg]
				const double magnet_longitude // �����̌o�x [deg]
			);
			axis_info( const axis_info& r );

			// ���̍Đݒ�
			// �����̈ܓx�o�x�́A�o�Ɏq���[�����g�����̕����ƂȂ锼�����̓_���L�q����B
			// �Ⴆ�΁A�n���̏ꍇ�͓씼����̎��ɂ�n���Ȃ���΂Ȃ�Ȃ��B
			void setAxis(
				const double magnet_latitude,
				const double magnet_longitude // �����̌o�x [deg]
			);

			// ���̖k�������̈ʒu���A�ܓx�E�o�x�n�ŕ\�����O�����������W�ŕԂ��B
			const vector& getGeometricRotationalAxis() const
			{ return m_rotAxis; }

			const vector& getGeometricMagneticalAxis() const
			{ return m_magAxis; }


		private:
			vector m_rotAxis;
			vector m_magAxis;
		};

	protected:
		// �R���X�g���N�^�E�f�X�g���N�^
		basic_planet(
			const double         radius,        // ���̘f���̔��a���w�肷��B
			const double            VDM,        // ���̘f���̉��z���C���[�����g���w�肷��B
			const axis_info&       axis,        // �f���̎��]�E��������n���B
			basic_magnet_model&     mag,        // ���ꃂ�f���̃C���X�^���X���w�肷��B
			basic_plasma_model&    plsm         // �v���Y�}���f���̃C���X�^���X���w�肷��B
		);

	private:
		void create();

	public:
		virtual ~basic_planet(){};

		
	public:

		// ���f���ւ̃A�N�Z�X ----------------------------------------
		basic_magnet_model& getMagnet() { return m_magnet; }
		const basic_magnet_model& getMagnet() const { return m_magnet; }

		basic_plasma_model& getPlasma() { return m_plasma; }
		const basic_plasma_model& getPlasma() const { return m_plasma; }


		// ���W�ϊ��n ------------------------------------------------
		// ��]�s�� --------------------------------------------------
		virtual matrix getGEI2GEO() const;
		virtual matrix getGEI2GSE() const;
		virtual matrix getGSE2GSM() const;
		virtual matrix getGSM2SM () const;


		// �n�����ǂ����̔��� ----------------------------------------
		double getRadius () const
		{ return m_radius; }

		bool isUnderSoil ( const vector& p ) const
		{ return restToSoil(p) <= 0.0; }

		virtual double restToSoil( const vector& p ) const
		{ return norm_2(p) - m_radius; }


		// �f������
		void setMagneticalAxis(
			const double magnet_latitude, // �����̈ܓx [deg]
			const double magnet_longitude // �����̌o�x [deg]
		);
		double getVirtualDipoleMagnet() const
		{ return m_vdm; }


		// �f���ɑ΂���ʒu���W�n ------------------------------------
		
		// �w��ʒu���玥�͐����g���[�X���A�t�b�g�v�����g���W��Ԃ��B
		// ���̃��\�b�h�́A���ꃂ�f���� getFootPrint()���Ăяo��
		// ���̌��ʂ�Ԃ��B
		vector getFootPrint(
			const vector& source_ptr,
			double      trace_factor
		) const;

		// �w��ʒu���玥�͐����g���[�X���A���C�ԓ���(Z=0)��̓_��Ԃ��B
		// ���̃��\�b�h�́A���ꃂ�f���� getEquatorPrint()���Ăяo��
		// ���̌��ʂ�Ԃ��B
		vector getEquatorPrint(
			const vector& source_ptr,
			double      trace_factor
		) const;
		
		// getFootPrint()�̌��ʂ���A�t�b�g�v�����g�ܓx(FLAT)��Ԃ��B
		// ���ʂƂ��āA�s�ώ��C�ܓx(ILAT)��Ԃ��B
		double getFLAT(
			const vector&  point,
			double  trace_factor
		) const;

		// getEquatorPrint()�̌��ʂ���AL�l�𓱏o���A
		// �o�Ɏq��������肵����Ŏ��C�ܓx��Ԃ��B
		// ���ʂƂ��āA�s�ώ��C�ܓx(ILAT)��Ԃ��B
		double getEqLAT(
			const vector&  point,
			double  trace_factor
		) const;

		// �o�Ɏq��������肵����ŁA���C�o�x�� 0����24�̐��l�ɒ����Ă������B
		// ���ʂƂ��Ď��C���[�J������(MLT)���������B
		double getMLT(
			const vector& point
		) const;
		
		// ���C�ܓx�E�o�x�ƍ��x��񂩂玥�͐����g���[�X���A�ړI�n�_�̍��W��Ԃ��B
		// ���̃��\�b�h�ł́A�n���ꂽ�ܓx�E�o�x�𖞂����n�\�ʂ̈ʒu����
		// ���ꃂ�f���ɉ����������Ɏ��͐����g���[�X���A�w��̍��x�ɂ��ǂ蒅�����_��
		// ���[�N���b�h��ԍ��W�n�ŕԂ��B
		vector getLocation(
			const double MLAT,      // MLAT�l���w��[degree] -90 <= MLAT <= 90
			const double MLT,       // MLT�l���w�� [h]        0 <= MLT  <  24
			const double altitude,  // �f���\�ʂ���̍��x���w��
			const double trace_factor = 1e3, // ����g���[�X�̐��x
			std::vector<vector>* const out_trace_line = NULL// �g���[�X�o�H���i�[
		) const;


	protected:
		// ���z�����A���]���A���� -------------------------------------
		// �����̃��\�b�h�Ɍ���AGSE���W�n�ŕԂ��B
		const vector  getRotationalAxisInGSE() const;
		const vector  getMagneticalAxisInGSE() const;


	protected:
		basic_magnet_model& m_magnet;
		basic_plasma_model& m_plasma;

	private:
		const double m_radius;        // �f���̔��a
		const double m_vdm;           // �f���̉��z���C���[�����g

		axis_info m_axisInfo; // ���̈ʒu

#	ifndef NDEBUG
		// �e�X�g�p�̃��\�b�h
		void test() const;
#	endif
	};
}

#endif//RTC_RAYTRACE_BASIC_PLANET_H

