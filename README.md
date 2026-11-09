# Goppa符号復号の実験コード

このリポジトリは、Goppa符号を中心とした符号理論・復号アルゴリズムのC言語実験コードの集まりです。
単一の完成品アプリケーションではなく、複数の独立した試作コードをまとめた研究向けのコードベースです。

主なテーマは、有限体上の演算、多項式の操作、そして誤り訂正符号の復号手順の実装です。
実験の対象としては、Berlekamp–Massey法、Euclid法、Patterson法、Peterson法などが含まれています。

## これは何をしているか

このコード群の目的は、Goppa符号に対する復号法を理論と実装の両面から検証し、各手法を比較することです。
一般に、Goppa符号の復号には次のような要素が関わります。

- シンドローム計算
- 誤り位置候補の決定
- 誤り位置多項式の生成
- 誤り値の解決
- 有限体や多項式環の演算

このリポジトリでは、それらを複数の別々の実装で試しています。

## ディレクトリの役割

| ディレクトリ | 概要 |
| --- | --- |
| `bm_gf2m/` | GF(2^m)上でのBerlekamp–Massey法の実験。 |
| `daikoku/` | Euclid法・Patterson法に近い復号の試作。 |
| `gf2m-peter/` | Peterson風の符号復号器・試験コード。 |
| `marimo/` | 拡大体上のパターソン復号アルゴリズム。 |
| `marimo/bunka2/`, `marimo/gaal/` | `marimo` の別バージョンや派生実装。 |

ルート直下にも、これらに関連する試作コードやヘッダがあります。ディレクトリ内のREADMEとMakefileを参照してください。符号パラメータは主に `global.h` または `global-p.h` で定義されています。

## CMakeでのビルド

リポジトリのルートでCMakeを実行すると、登録済みの実験プログラムをまとめてビルドできます。

```sh
cmake -S . -B build
cmake --build build --parallel
```

現在登録されているターゲットは、各コード群の代表的な実行ファイルを個別にビルドするためのものです。リポジトリ内のすべての `.c` ファイルを一括でコンパイルするプロジェクトではありません。

| CMakeターゲット | ソース |
| --- | --- |
| `tokyo` | `ymo.c` |
| `bm_gf2m_lyon` | `bm_gf2m/fat_old_sun.c` |
| `daikoku_xeuclide` | `daikoku/Xeuclide.c` |
| `gf2m_peter_ptar` | `gf2m-peter/origin.c` |
| `marimo_vmo` | `marimo/vmo.c` |

実行したいターゲットを個別に起動します。

```sh
./build/marimo_vmo
```

`daikoku_xeuclide` と `marimo_vmo` はAVX2向けビルドを含みます。
そのため、実行環境のCPUがAVX2対応であることが必要です。

## Niederreiter型デモ

`marimo/niederreiter_demo.c` はBM/Patterson実装を使った往復確認用デモです。
リポジトリのルートから次のようにビルド・実行します。

```sh
make -C marimo niederreiter_demo
./marimo/niederreiter_demo [誤り重み] [raw]
```

誤り重みは `1` から `K` の範囲で指定します。省略時は `K`（現在の設定では16）です。
通常はランダムな行変換と列置換を行います。第2引数に `raw` を指定すると両方を無効にし、
スクランブル前の行列で比較します。`raw` を使う場合も誤り重みを第1引数に置いてください。

```sh
./marimo/niederreiter_demo 1       # 誤り重み1、スクランブルあり
./marimo/niederreiter_demo 1 raw   # 誤り重み1、スクランブルなし
./marimo/niederreiter_demo 16 raw  # 既定の誤り重みでスクランブルなし
```

`raw` だけを第1引数にすると誤り重みとして解釈され、引数エラーになります。
このデモの詳細と出力の読み方は `marimo/README.md` を参照してください。

## 注意事項

- このコードは研究・検証・学習用途の実験コードです。
- 完成した実用ライブラリではなく、理論検証のための試作コードです。
- 実装の状態や前提条件はディレクトリごとに異なります。
- ルートの従来Makefileは古い依存関係を持つため、CMakeを使う方式がより安定しています。

## まとめ

このリポジトリは、「Goppa符号の復号法を理論から実装まで試すための集合的な実験コード群」です。
復号法の違いを比較しながら、有限体演算や多項式手法を調べるためのコードとして設計されています。

このコードは研究・学習・アルゴリズム検証向けの試作です。実装の完成度や暗号学的安全性は保証されていません。実用の暗号システムや機密データの保護には使用しないでください。

## 今回の肝とでも言えるコード。
**BM用のシンドロームは、Goppa多項式の係数行列をかけることでパターソン用のシンドロームへと変換できる。**

```c
// synd() の vec.x に格納された通常順シンドロームを、g の係数から作る
// 上三角 Toeplitz 行列で Patterson 用シンドロームに変換する。
static q_res q_bm_to_patterson_syndrome(const short *syndrome)
{
    static short toeplitz[K][K];
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            toeplitz[i][j] = (j >= i) ? q_g[K + i - j] : 0;

    q_res s = q_zero();
    for (int i = 0; i < K; i++)
    {
        short acc = 0;
        for (int j = i; j < K; j++)
            acc = q_add(acc, q_mul(toeplitz[i][j], syndrome[K - 1 - j]));
        s.c[i] = q_neg(acc);
    }
    return s;
}
```
