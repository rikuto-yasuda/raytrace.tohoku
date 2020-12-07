
RAYTRACE='./testing'


## パラメータを設定 ##############################
## 発生源を指定 ##
# euclid指定の場合は、惑星半径単位で。
# polar指定の場合は、MKS単位で。
#COORD="polar"   # (euclid|polar)
COORD="euclid" # (euclid|polar)
#SX=60           # (source.x|MLAT)
#SX=70           # (source.x|MLAT)
#SY=23           # (source.y|MLT)
#SZ=2e6         # (source.z|altitude)
#SZ=1.1e8        # (source.z|altitude)

SX=3
SY=4
SZ=9


## モデルを選択 ##
PLASMA="test_simple"    #(null|test_null|simple|test_simple|sato|nsumei|devine_garrett)
MAGNET="test_simple"         #(null|test_null|simple|test_simple|igrf|igrf4|vip4)
PLANET="jupiter"        #(earth(?)|jupiter|benchmark)

## 宇宙の時刻を指定 ##
DATE="2000/1/1"  # year/month/day
TIME="0:0.0"     # hour:minutes.sec

## 波動特性を設定 ##
#FREQ=0          # 周波数[Hz]
FREQ=30e6       # 周波数[Hz]
MODE="LO"       # 波動モード(LO|RX)
#RAY_L=1e6       # トレースする最大の光路長
RAY_L=7.1e8     # トレースする最大の光路長
PITCH=90        # 磁場に対するピッチ角
SEGMENT=30      # 出力する光路上の点の数
MAX_STEP=5000   # トレース・ステップの最大数
STEP_LENGTH=1e8 # １ステップで進む最大の光路長
#STEP_LENGTH=1e7 # １ステップで進む最大の光路長
PRECISION="3.74e-5"  # １ステップ間のベクトル誤差の許容率
#PRECISION="3.74e-6"  # １ステップ間のベクトル誤差の許容率
TIME_RANGE="1:1e-6"  # １ステップ間の時間分解能レンジ
#TIME_RANGE="1:1e-6"  # １ステップ間の時間分解能レンジ

## plasma cavity ##
# --cavity [fp/fc]/[ilat]:[ilat range]/[mlt]:[mlt range]/[height upper]:[height bottom]
# すべてMKSA単位系
CAVITY_LIST=(                      \
  '--cavity 0.03/70:3/0:1/3e7:5e6' \
) # cavityの数だけオプションを指定

## 出力ファイル名を指定する。
OUTPUT="ray-P${PLASMA}-M${MAGNET}-${PLANET}-${MODE}"
LOG="${0}.log"

## 終了時にメールを送る ##
MAIL_TO="" # 終了メールを送る先を指定する。
MAIL_FROM="${USER}" # ここでFromを明示することができる。
MAIL_SUBJECT="[mkray] ${OUTPUT} was completed."


##################################################

send_mail()
{
	if [ ${MAIL_TO} ]; then
		echo "${MAIL_SUBJECT}" | mail "${MAIL_TO}" -s "${MAIL_SUBJECT}" -- -f ${MAIL_FROM}
	fi
}

## main ##########################################
	echo "BEGIN (${OUTPUT}) at " `date` >> ${LOG}

	$RAYTRACE \
	  --plot ray                    \
	  --verbose                     \
	  --source-coord ${COORD}       \
	  --source-x     ${SX}          \
	  --source-y     ${SY}          \
	  --source-z     ${SZ}          \
	  --plasma-model ${PLASMA}      \
	  --magnet-model ${MAGNET}      \
	  --freq         ${FREQ}        \
	  --ray-mode     ${MODE}        \
	  --ray-length   ${RAY_L}       \
	  --step-count   ${MAX_STEP}    \
	  --step-length  ${STEP_LENGTH} \
	  --time-range   ${TIME_RANGE}  \
	  --precision    ${PRECISION}   \
	  --pitch        ${PITCH}       \
	  --ray-path-segment ${SEGMENT} \
	  --planet       ${PLANET}      \
	  ${CAVITY_LIST}                \
	  2>&1                          \
	  1> ${OUTPUT}                  \
	  | tee ${LOG}

#使っていないオプション
#	  --back-trace                  \
#	  --without-plot-startptr       \
#	  --parallel                    \


	echo "END (${OUTPUT}) at " `date` >> ${LOG}
	send_mail
