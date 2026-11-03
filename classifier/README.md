# classifier/ — SILVA 学習済み分類器の置き場

🔴 **`.qza` と一緒に `README`・`.doi` を必ず置くこと。** パイプラインは実行時にこの 2 つを読み、
分類器の素性（SILVA 版・QIIME 2 版・領域・habitat・DOI）をログに残す（config には書かない）。

`.qza` 本体は **Git 管理・バックアップの対象外** (261 MB・公式から再取得可能)。
`.md5`・`.doi`・SILVA 同梱 `README` は**本体と必ず一緒に置く** (再現性と引用のため)。

## レイアウト

```
classifier/
├── README.md              ← この文書
├── uniform/               ★ 既定 (属まで・classic)
│   └── ..._uniq_genus_uniform_classifier.qza   261 MB  ✅
├── average/               ○ 任意 (種まで・14 habitat の平均)       904 MB  ✅
├── human-stool/           ○ 任意 (種まで・ヒト糞便)                903 MB  ✅
└── animal-distal-gut/     ○ 任意 (種まで・動物糞便=マウスに最も近い) 903 MB  ✅

**260913: 4 種すべて取得・md5 照合済み（計 2.9 GB）。** 各フォルダに `.qza` / `.md5` / `.doi` / SILVA 同梱 `README`。
**既定は `uniform`**（決定の経緯は `../README.md` の B-5b）。weighted は config 差し替えだけで使える。
```

各 weighted フォルダは `.md5` / `.doi` / `README` のみ取得済み。**照合すべき md5**:

| habitat | md5 |
| --- | --- |
| average | `49129868a1d12c413879f88b28e05429` |
| human-stool | `76066593508e0073873455246a574e8a` |
| animal-distal-gut | `6a88eb76b15a6af3aea12913d4376083` |

### 動作確認 (260912)

```
$ qiime tools peek <uniform の .qza>       # rachis-qiime2-2026.7
UUID:        b7986601-1851-4fa2-ace5-98bacf660fa3
Type:        TaxonomicClassifier
Data format: TaxonomicClassiferTemporaryPickleDirFmt
```
→ **QIIME 2 2026.7 で正しく読める artifact であることを確認済み。**

取得元: https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/full-length/uniform/

### 同フォルダの 3 ファイルの区別 (260912 実物照合)

| ファイル | サイズ | 判断 |
| --- | --- | --- |
| `..._rank-propagation_uniq_genus_uniform_classifier.qza` | 261 MB | ✅ **本体。これを使う** |
| `SILVA_144_SSURef_NR99_uniform_classifier_full-length.qza` | 261 MB | ❌ **同一物** (md5・DOI とも一致)。旧命名の別名と推測。落とさない |
| `rescript-evaluation_full-length.tar.gz` | 21 MB | 🔵 任意・推奨。**分類器の精度レポート** |

**ファイル名の読み方**: `SSURef_NR99` = SSU (16S/18S) の高品質版を 99% で重複除去 /
`rank-propagation` = 空欄の階級を上位名で穴埋めして階級を揃えてある /
`uniq` = 同一配列を除去 / **`genus`** = 属まで (weighted は `species`) / `uniform` = 出現確率を一様と仮定。

**rescript-evaluation の中身** (SILVA 同梱 README より・3 つの .qzv):
`fit-classifier-evaluation` (学習の評価) / `classifier-evaluation` (**階級ごとの precision・recall・F 値**) /
`slv-taxonomy-evaluation` (分類体系そのものの統計)。**解析には使わないが、査読で
「なぜ属までなのか」「分類の信頼度は」と問われたときに SILVA 自身の精度指標を引ける。**
md5 `7e6ffbd2a566e4a3df48c7686f831535` / DOI `.../full-length/uniform/evaluation`
※ 各 qzv が具体的に何を描くかは RESCRIPt の一般的挙動からの推定。開けば確定する。

## 照合すべき md5

```
e16b9bc78dafa2f344563e570bbf25ba
```

```bash
cd <パイプラインのフォルダ>/classifier/uniform && md5 -r *.qza && cat *.qza.md5
```

⚠️ **パイプラインは md5 照合を行わない**（260915 判断で廃止・config からも `md5` キーを削除）。
**ダウンロード直後に、この手順で人が 1 度だけ確認する。**

## DOI (論文の Methods に書く)

```
10.82364/144/SSU/Ref-NR99/QIIME2/2026.7/full-length/uniform
```

## 引用文献 (SILVA 同梱 README より)

- **SILVA**: Quast C et al. (2013) *Nucleic Acids Res* 41(D1):D590-D596
- **QIIME 2**: Bolyen E et al. (2019) *Nat Biotechnol* 37:852-857
- (weighted を使う場合のみ) **q2-clawback**: Kaehler BD et al. (2019) *Nat Commun* 10:4643 /
  Bokulich NA et al. (2018) *Microbiome* 6:90

## GPG 署名の検証 (md5 より強い確認・任意)

```bash
qiime tools signature-verify --input-path <artifact.qza> --name 'SILVA-Pipeline-Official'
```
公開鍵 `SILVA_PUBLIC_KEY.asc` を `gpg --import` してから実行する。

## ⚠️ SILVA 同梱 README の落とし穴

"How to Use" のコマンド例にある `silva-138-99-515-806-nb-classifier.qza` は
**旧版の名前の消し忘れ**。実際には本フォルダの 144 の `.qza` を指定する。

## 任意 (weighted・種まで分類・experimental)

`.../SSU/full-length/weighted/{average,human-stool,animal-distal-gut,...}/`
⚠️ **habitat が検体と違うと誤った結果になる** (SILVA 公式の注意)。
マウス糞便なら `animal-distal-gut` が最も近いが、**既定は uniform のまま**とし、
weighted は併用検証で使う。

## 取らないもの

`/current-release/QIIME2/2026.7/taxonomic-weights/` は **分類器を自分で焼くための材料**
(habitat 別の出現頻度データ = q2-clawback の入力)。上の weighted 分類器は SILVA が
この重みで焼いた完成品なので、材料側は不要。
