# 16S 解析パイプライン 使い方の手引き

MiSeq の 16S rRNA 遺伝子 (V1-V2 領域・27Fmod / 338R) の paired-end FASTQ から、QIIME 2 で ASV の表・多様性・分類を作り、群の比較の図 (#1〜#12) まで描くパイプラインです。1 つのコマンドで 01〜10 の手順を通して実行します。

- 版: `VERSION` に書いてあります。変更の記録は `CHANGELOG.md` です。
- 設定の書き方の詳細は、雛形 `config.example.yaml` のコメントにも書いてあります。
- このパイプラインは、内部標準 (IS) の配列とネガティブコントロールの除去を行いません。
  IS を入れた検体では IS の ASV が表に残り、相対存在量や多様性の値に入ります。
  ネガティブコントロールは、ほかの検体と同じに扱われます。解析に入れないときは、config の `exclude` に書きます (§4)。

---

## 0. 必要なもの

| もの | 内容 |
| --- | --- |
| QIIME 2 2026.7 (amplicon) の conda 環境 | 公式の手順で作った env。R と R のパッケージ (dada2・ANCOMBC・phyloseq)、scikit-bio などはこの env に含まれています。env の名前は何でもかまいません |
| このパイプラインのフォルダ一式 | `run_16s.py`・`lib/`・`config.example.yaml`・`classifier/`・この手引きなど |
| SILVA 144 の分類器 | `classifier/` に同梱 (§2) |
| 空き容量 | 分類器 2.9 GB のほか、実験 1 つにつき FASTQ と出力で数 GB (126 検体のランで FASTQ 2.1 GB・出力 2.3 GB) |
| 時間 | 126 検体・1 ランで約 1 時間 (Apple M2・メモリ 8 GB での実測)。DADA2 と 09 の分類に時間がかかります |

> **[未記入] OS ごとの環境の作り方**: 使う OS (Windows の WSL・Linux・macOS) が決まったら、QIIME 2 の公式の手順へのリンクと、確かめ方をここに書きます。

**env の確かめ方** (env に入った状態で):

```bash
qiime info                       # QIIME 2 version: 2026.7.x と出ること
Rscript -e 'suppressMessages({library(ANCOMBC); library(phyloseq)})' && echo OK   # 最後に OK と出ること
```

パイプラインも起動するたびに、QIIME 2 の版・R の場所 (`R_HOME`)・ANCOMBC と phyloseq を確かめ、足りなければ理由を出して止まります。

---

## 1. 置き方

パイプラインのフォルダと、実験ごとのフォルダ (実験フォルダ) を分けて置きます。実験フォルダの場所と名前は自由です (コードは名前を見ません)。名前の決まりがあるのは、ランのフォルダと検体情報 CSV だけです。

```
<実験フォルダ>/
├── input/
│   ├── 16S/
│   │   ├── config.yaml                  雛形 config.example.yaml を写して書く (§4)
│   │   └── runs/
│   │       └── <YYMMDD>_M<4 桁>/          ランのフォルダ (シーケンスした日 + ラン番号。例 250101_M0001)
│   │           ├── fastq/                 FASTQ だけを入れる (§3.1)
│   │           ├── SampleSheet.txt        そのランの検体表 (§3.2)
│   │           └── rename.tsv             古い名前を付け直したときだけ (§3.4)
│   └── sample_information/
│       └── <YYMMDD>_<P 番号 4 桁>_<名前>.csv   検体情報 CSV (§3.3)
└── output/16S/                            パイプラインが作る (§6)
```

- ランの FASTQ は、そのランを最初に使った実験フォルダに置きます。ほかの実験でそのランを使うときは、コピーせず config の `runs:` にパスを書きます。
- config に書くパスは、すべて `~` で始まるパスか絶対パスにします (相対パスは止まります)。フォルダを移したら、config のパスも書き直します。

---

## 2. 分類器 (`classifier/`)

SILVA 144 SSURef NR99 の学習済み分類器 (QIIME 2 2026.7 用・full-length) を 4 種同梱しています。

| フォルダ | 内容 | 使い方 |
| --- | --- | --- |
| `uniform/` | 属まで・出現確率を一様と仮定 (261 MB) | **既定** |
| `average/`・`human-stool/`・`animal-distal-gut/` | 種まで・habitat ごとの出現確率で重み付け (各 903 MB) | 任意。SILVA は「experimental・検体と habitat が違うと誤った結果になる」と注意しています |

- 各フォルダの `.qza` の隣にある `README` (SILVA 同梱) と `.qza.doi` は消さないでください。パイプラインはこの 2 つを読んで、分類器の素性 (SILVA・QIIME 2 の版・領域・DOI) をログに残します。
- config の `classifier.path` に、使う `.qza` の場所を書きます (§4)。
- 取り直すときの入手先: `https://www.arb-silva.de/current-release/QIIME2/2026.7/SSU/full-length/` (uniform の md5 は `e16b9bc78dafa2f344563e570bbf25ba`)。パイプラインは md5 を確かめないので、ダウンロードしたら一度だけ手で確かめてください。

---

## 3. 入力の作り方

### 3.1 FASTQ (`runs/<ラン>/fastq/`)

- 名前: `<検体の ID>_S<番号>_L001_R1_001.fastq.gz` と `…_R2_001.fastq.gz` の組 (例 `M0001P0001N0001_S1_L001_R1_001.fastq.gz`)。
- 検体の ID は `M<ラン番号 4 桁>P<プロジェクト番号 4 桁>N<通し番号 4 桁>` の形だけです。
  - `M` = ランの番号 (ランのフォルダ名の M と同じ)
  - `P` = 検体情報 CSV の `Project No.`
  - `N` = 検体情報 CSV の `No.`
- タグ分けした検体は、`S` の欄を `Ssplit` にします (例 `M0001P0001N0001_Ssplit_L001_R1_001.fastq.gz`)。
- フォルダには、自分の解析で使う検体の FASTQ だけを入れます。次のものがあると止まります。
  - FASTQ 以外のファイル (Mac の `.DS_Store`・`._…` などの隠しファイルを含む)
  - `Undetermined_…` の FASTQ (MiSeq が出す、どの検体にも振り分けられなかった読み)
  - R1 と R2 がそろっていない検体、ほかのランの番号の FASTQ
- 新しいランは、供給元にこの名前で納品してもらいます。古い形の名前 (`SO001` など) は、手元で付け直してから使います (§3.4)。

### 3.2 SampleSheet.txt (`runs/<ラン>/SampleSheet.txt`)

- タブ区切り・UTF-8 のテキスト。1 列目の列名は `#SampleID` です。
- 2 列目以降 (例 `BarcodeSequence`・`LinkerPrimerSequence`・`ReversePrimer`・`Description`) は、そのまま検体の表に写されます。
- ラン全体の表のままでかまいません (相乗りしたほかの方の検体の行があってもよい)。FASTQ のある検体の行だけを使います。
- 止まる場面: FASTQ の ID が `#SampleID` に無い、同じ `#SampleID` が 2 行以上ある (使わない行でも)、UTF-8 で読めない。

### 3.3 検体情報 CSV (`input/sample_information/`)

- ファイル名: `<YYMMDD>_<P 番号 4 桁>_<名前>.csv` (例 `260709_0503_Sample information document.csv`)。
- 文字コード: UTF-8 (Excel なら「CSV UTF-8」で保存)。UTF-8 のエラーが出たときは §8 を見てください。
- 形 (行の位置は固定。ひな形 `templates/sample_information_template.csv`):

| 行 | 中身 |
| --- | --- |
| 1 | `Sample information` |
| 2 | `Project No.,<番号>` (ファイル名の P 番号と同じ数。例 503) |
| 3〜5 | プロジェクト名・日付・担当者など (パイプラインは読まない) |
| 6 | 空行 |
| 7 | 列名。`No.`・`Sample_name`・`Group_name`・`Subgroup_name` が必須 (ほかの列があってもよい) |
| 8〜 | 1 行 1 検体 |

- `No.` = 整数 (検体の ID の `N`)。重複は不可。
- `Sample_name` = 検体名。図や表の検体の名前 (04 以降の ID) になります。全 CSV を通して重複は不可。
  数字以外の文字を 1 文字以上入れてください (`727` は不可、`_727`・`M727` は可)。
  QIIME 2 は数だけの列を数量として扱うので、検体名に使えません。
  数として読める名前が 1 つでもあると、01 で止まります。
  図 #7 (積み上げ棒・検体ごと) では、名前が `_番号` で終わると、群の中を番号の順に並べます。
  それ以外は文字の順です (`M10` が `M9` より前)。桁数をそろえれば番号の順と同じになります。
- `Group_name` = 群 (図の色・検定の単位)、`Subgroup_name` = サブグループ (基質など。図の記号)。
- 使える文字: 半角の英数字と `. - _ : |` だけ (空白・`+`・`%`・`&`・全角は不可)。36 文字以内で、`#` で始めない。`id`・`sampleid` などの予約語は使えません。
  例: `0.1% Diet A + B` → `DietA-B`
- CSV に載っていて FASTQ の無い検体 (MiSeq に出していない検体) は、config の `absent_ok` に書きます (§4)。書いていないと、入れ忘れとみなして止まります。
- ほかの列 (採取日・量など) は、図の版のフォルダの `00_info/sample_information.csv` に写されます。

### 3.4 古い名前の FASTQ を付け直すとき (`rename.tsv`)

- FASTQ の名前の ID を新しい形に付け直し、SampleSheet.txt の `#SampleID` も同じ名前に書き換えます。
- 付け直しの対応表 `rename.tsv` (タブ区切り。1 行目 `old	new`、以下 `旧ファイル名	新ファイル名`) をランのフォルダに残します。パイプラインは読みませんが、元の名前を辿るための記録です。

---

## 4. config の書き方 (`input/16S/config.yaml`)

雛形 `config.example.yaml` を写して、値だけを書き換えます (コメントは残してください)。**すべての項目が必須**で、使わないものは `none` と書きます。書き漏れ・書き間違いは、01 の前に止まって知らせます。記入済みの例: `templates/config_example_filled.yaml`。

| 項目 | 書くこと |
| --- | --- |
| `threads` | 並列数。0 = 全コア |
| `runs:` | 使うランのリスト。各ランに `path` (ランのフォルダ)・`quality_control`・`exclude` |
| `quality_control` | `trim_left_f`・`trim_left_r` = 先頭から切る長さ (プライマーの長さ。27Fmod 20・338R 19。タグ分けのランは 23・22)、`trunc_len_f`・`trunc_len_r` = 品質が落ちる手前で切る長さ (§5 の `--stage 1` で決める) |
| `exclude` | そのランで解析に入れない検体の ID のリスト (例 `[M0001P0001N0015]`)。無ければ none |
| `parent:`・`select:` | 絞り込みの形 (§7) のときだけ。ふつうは両方 none |
| `sample_info:` | 検体情報 CSV のリスト。各 CSV に `path` と `absent_ok` (FASTQ が無くてよい検体の No.。例 `[N0127, N0128]`。無ければ none) |
| `diversity` | `sampling_depth_floor` = 希薄化の深さの下限 (研究室の慣習は 10000)、`random_seed` = 希薄化の乱数の種 (0 のままでよい) |
| `classifier` | `path` = 分類器の `.qza` の場所、`n_jobs` = 分類の並列数 (1 つにつき分類器 1 つ分のメモリを使う。メモリ 8 GB なら 2。初回のあと `09_taxonomy/classify_peak_rss_mb.txt` を見て調整) |
| `taxonomy` | `confidence` (0.7 のままでよい)・`collapse_ranks` (表にする階級) |
| `grouping` | 群の並び・色・サブグループの記号 (下の注意) |
| `comparison`・`stats` | 対照群・検定の方法 (tukey・dunnett・welch)・検定に要る最小の n |
| `plotting`・`pipeline.run` | 図に載せる数、どの図を描くか |

**grouping の注意**
- `group_order` には、データにある群をすべて書きます (無い群・余分な群があると止まります)。
- `use_subgroup: true` にすると、群が「群 サブグループ」の組み合わせになり、検定も変わります。対照群も `"Control SubA"` のように書きます。
- `show_group_subgroup: true` は、色は群のまま、サブグループの記号を重ねるだけです (検定は変わりません)。
- 色は、全群を個別に書く (`"#228833"` の形) か、配色名 (`muted`・`bright`・`dark`・`smooth_rainbow` など) か、none です。

**複数のランを使うとき**: `runs:` にランを並べます。ランごとに DADA2 を回したあと、1 つの表につなげます。trim_left は、切ったあとの開始位置がランどうしで同じになるように書いてください (タグ分けのランは +3)。同じでないと、同じ菌が別の ASV に分かれ、ランの差が群の差のように見えます。パイプラインはこれを確かめません。

---

## 5. 実行

```bash
conda activate <QIIME 2 2026.7 の env>
python <パイプラインのフォルダ>/run_16s.py <実験フォルダ>
```

- `<パイプラインのフォルダ>` は、zip を展開した場所です。古いコピーが残っていると、古い版が動きます。
  実行した版と場所は、ログの冒頭の `パイプライン :` の行で確かめられます。
- 実験フォルダを省くと、起動してから尋ねます。
- 最初は `--stage 1` を付けて読み込み (02) だけで止め、`output/16S/Qiime2/HTML_02_demux_<M>.qzv/index.html` の Interactive Quality Plot を見て `trunc_len` を決めるのがおすすめです。ルーチンの MiSeq なら、雛形の値 (280・210) で問題ないことが多いです。
- 途中で止まっても、同じコマンドでやり直せます。出力ができている手順は飛ばして、続きから進みます。
- 実行中は、パソコンをスリープさせないでください (1 時間ほどかかります)。

### 5.1 端末で打つときの注意

- 引用符はまっすぐな半角の `"`・`'` を使います。`“ ”`・`‘ ’` は使えません (メールや文書から写すと変わることがあります)。
- 空白を含むファイル名は、Tab キーで補うと確実です。
- 絶対パスは `/` から書きます (`home/…` は不可、`/home/…` は可)。
- ssh でつないで実行するときは、`tmux` の中で実行します。接続が切れても実行が続きます。
  `conda activate` も `tmux` の中でやり直します。

  ```bash
  tmux new -s qiime
  ```

  つなぎ直したあとは `tmux attach -t qiime` で戻れます。
- 画面の無い計算機では、`.qzv` ファイルを手元のパソコンに写し、https://view.qiime2.org で開きます。
  写すのは番号のフォルダの中の `.qzv` です (例 `09_taxonomy/taxa-bar-plots.qzv`)。
  `HTML_…qzv` はフォルダなので、この方法では開けません。

---

## 6. 出力 (`<実験フォルダ>/output/16S/`)

```
output/16S/
├── Qiime2/          01〜09 の QIIME 2 の出力と、HTML_<手順>_<名前>.qzv/ (ブラウザで index.html を開く)
│   ├── 01_metadata/       検体の表 (sample-metadata.txt・sample-metadata_cn.txt・sample_information.csv)
│   ├── 02_import/<M>/     ランごとの読み込み
│   ├── 03_dada2/<M>/      ランごとの DADA2 (03_merge/ = ランが 2 つ以上のときの結合)
│   ├── 04_relabel/〜09_taxonomy/
├── figure/ver<N>/   図。00_info/ (凡例・設定・名前の対応表・検体情報の表) と、図の番号ごとのフォルダ (PDF と csv/)
└── logs/            実行ごとのログ ({日時}_16S.log) と env の一覧 ({日時}_env.txt)
```

**図の版 (`ver<N>`)**: 図の設定 (grouping・comparison・stats・plotting・絞り込み) か、上流の条件が前と違えば新しい版のフォルダができます。同じなら同じ版に描き足します (図のフォルダがある番号は飛ばします)。描き直したい図は、そのフォルダをゴミ箱へ移してから実行します。

**図の一覧**

| # | 図 | 見るもの |
| --- | --- | --- |
| 1 | 最終リード数と sampling depth | 検体ごとのキメラ除去後のリード数。白抜き = 希薄化の深さに届かず外れた検体 |
| 2 | レアファクションカーブ | 読みを増やしたときに多様性が頭打ちになるか (深さが足りているか) |
| 3 | α多様性 | 群ごとの多様性 (observed features・Shannon・Faith PD・evenness) と検定 |
| 4 | PCoA | 群集の似かよりの散布図 (4 つの距離。PCo1×PCo2 と PCo1〜4 の組み合わせ) |
| 5 | PERMANOVA / PERMDISP | 群の組ごとの、中心の違い (PERMANOVA) とばらつきの違い (PERMDISP) |
| 6 | デンドログラム | 検体ごと・群平均の階層クラスタリング |
| 7 | 積み上げ棒 | 群平均・検体ごとの組成 (上位のタクソン・残りは Others) |
| 8 | ヒートマップ | 上位のタクソンの群平均 (対照群との差の星) |
| 9 | ベン図 / UpSet | 群ごとに検出されたタクソンの重なり (ベン図は 2〜5 群のときだけ) |
| 10 | タクソン別の棒 | 1 タクソン 1 枚の相対存在量と検定 |
| 11 | ANCOM-BC | 対照群に対して増えた・減ったタクソン (ANCOM-BC2) |
| 12 | ANCOM-BC まとめ | #11 に、対照群での存在量の列を足した図 |

- 図と同じ値は、各フォルダの `csv/` にあります。
- 希薄化の深さに届かなかった検体は、#1 には出ますが、#3 以降の解析には入りません (`06_rarefaction/excluded_samples.tsv`)。

### 6.1 階級と番号の対応

同梱の SILVA 144 の分類器の番号です。k__ の階層が入るので、1 つずれます。

| `--p-level`・`level-N.csv` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 階級 | Domain | Kingdom | Phylum | Class | Order | Family | Genus | Species |

- 例: `level-6.csv` は科、`level-7.csv` は属です。`collapsed_table_genus.tsv` は `level-7.csv` と同じ階級です。
- `uniform` の分類器は、level 7 (Genus) までです。

---

## 7. 絞り込みの形 (一部の検体だけで図を見る)

解析を終えた実験 (親) の結果を使い、検体を選んで図 (10) だけを描き直せます。DADA2・希薄化・距離は親と同じ値のまま、PCoA・相対存在量・ANCOM-BC2 は選んだ検体で計算し直します。

1. 新しい実験フォルダを作り、`input/16S/config.yaml` だけを置きます。
2. config に `parent:` (親の実験フォルダ) と `select:` を書き、`runs`・`threads`・`sample_info`・`diversity`・`classifier`・`taxonomy` を none にします。

```yaml
parent: "~/<親の実験フォルダ>"
select:
  Group_name: all              # 使う群のリストか all
  Subgroup_name: [SubA]        # 使うサブグループのリストか all
  exclude: none                # さらに外す検体の ID のリストか none
```

- `grouping.group_order` には、選んだあとのデータにある群だけを書きます。
- 親が作り直されたら止まります (古い結果と混ざらないように)。
- 検体を最初から外して解析し直したい (DADA2 や希薄化も外したあとの検体で決めたい) ときは、絞り込みではなく、新しい実験フォルダで `runs:` の `exclude` に書いて 01 から実行します。

---

## 8. ログの読み方と、止まったとき

- ログは画面と `output/16S/logs/{日時}_16S.log` の両方に出ます。最後に「注記: N件」として、注意と注記がまとめて出ます。
- 印の意味:
  - `[エラー]` = 止まります。理由と、`※` で始まる直し方の案内が続きます。
  - `[注意]` = 図や結果に問題が出ているかもしれず、見直しが要ります (例 群名が見切れている・R の警告)。
  - `[注記]` = 想定どおりの処理の知らせで、対処は要りません (例 n が足りない群を検定から外した)。
- **解析の条件を変えるとき** (config の 01〜09 に関わる値・検体情報 CSV・FASTQ を変えたとき): 同じ実験フォルダで実行すると、前の条件で作った出力と食い違うので止まります (`[エラー] 既存の出力は、今の config と違う条件で作られています`)。**新しい実験フォルダを作って実行してください。** 前の条件のまま進めるなら、config を前の値に戻します。
- 図の設定 (grouping 以下) だけを変えたときは、同じ実験フォルダで実行してかまいません (新しい版の図ができます)。
- `Segmentation fault` や R の場所のエラーが出たら、`conda activate` をやり直してから実行します。
- `UTF-8 でない文字があります` と出たとき: 表示された行の `[ここ]` の文字を消すか打ち直し、UTF-8 (Excel の「CSV UTF-8」) で保存し直します。
  保存し直しても同じエラーが出るときは、別の文字コードの文字が混じったままです。
  多くの行が表示されるときは、ファイル全体が別の文字コード (Shift-JIS など) で保存されています。
  行番号はファイルの行で数えるので、セルの中に改行があると Excel の行番号とずれます。

---

## 9. 論文に書くとき

- 実行ログの冒頭に、パイプラインの版と実行した場所、QIIME 2 の env・作図に使ったライブラリの版 (scipy・scikit-bio・R ANCOMBC など)・分類器の素性 (SILVA の版・DOI) が残ります。Methods はこれを見て書きます。env の全パッケージの一覧は `logs/{日時}_env.txt` です。
- 分類器の DOI (uniform): `10.82364/144/SSU/Ref-NR99/QIIME2/2026.7/full-length/uniform`
- 引用文献:
  - QIIME 2: Bolyen E et al. (2019) *Nat Biotechnol* 37:852-857 (PMID 31341288)
  - SILVA: Quast C et al. (2013) *Nucleic Acids Res* 41(D1):D590-D596 (PMID 23193283)
  - DADA2: Callahan BJ et al. (2016) *Nat Methods* 13:581-583 (PMID 27214047)
  - ANCOM-BC2: Lin H & Peddada SD (2024) *Nat Methods* 21:83-91 (PMID 38158428)
  - (weighted の分類器を使う場合) Kaehler BD et al. (2019) *Nat Commun* 10:4643
- uniform を選んだときは「SILVA の注意書きに従い uniform の分類器を用いた」と書けます。
